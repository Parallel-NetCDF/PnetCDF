/*
 *  Copyright (C) 2026, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * Internal commux runtime for C callers.
 *
 * This follows the same split used by snapy: the internal runtime owns the
 * c10d TCPStore and commux::ProcessGroupUCX, while public bindings are thin
 * wrappers over this C ABI.
 */

#include <pnetcdf_commux.h>

#include <commux/process_group_ucx.hpp>
#include <torch/torch.h>
#include <torch/csrc/distributed/c10d/TCPStore.hpp>

#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <fcntl.h>
#include <cinttypes>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

namespace {

struct Runtime {
  at::intrusive_ptr<c10d::Store> store;
  c10::intrusive_ptr<commux::ProcessGroupUCX> pg;
  int rank = -1;
  int size = -1;
  int64_t next_file_handle = 1;
  std::map<int64_t, int> files;
};

std::mutex g_mutex;
std::unique_ptr<Runtime> g_runtime;
char g_last_error[1024] = "";

void set_error(const char* fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  std::vsnprintf(g_last_error, sizeof(g_last_error), fmt, ap);
  va_end(ap);
}

const char* getenv_required(const char* name) {
  const char* value = std::getenv(name);
  if (value == nullptr || value[0] == '\0') {
    set_error("missing required environment variable %s", name);
    return nullptr;
  }
  return value;
}

bool parse_tcp_init_method(const char* init_method, std::string& host,
                           int& port) {
  if (init_method == nullptr || std::strcmp(init_method, "env://") == 0) {
    const char* env_host = getenv_required("MASTER_ADDR");
    const char* env_port = std::getenv("PNC_COMMUX_MASTER_PORT");
    bool explicit_pnc_port = env_port != nullptr && env_port[0] != '\0';
    if (!explicit_pnc_port) env_port = getenv_required("MASTER_PORT");
    if (env_host == nullptr || env_port == nullptr) return false;
    host = env_host;
    port = std::atoi(env_port);
    if (!explicit_pnc_port) port += 101;
    return port > 0;
  }

  const char* prefix = "tcp://";
  const std::size_t prefix_len = std::strlen(prefix);
  if (std::strncmp(init_method, prefix, prefix_len) != 0) {
    set_error("unsupported commux init method '%s'", init_method);
    return false;
  }

  std::string endpoint(init_method + prefix_len);
  std::size_t colon = endpoint.rfind(':');
  if (colon == std::string::npos || colon + 1 >= endpoint.size()) {
    set_error("invalid commux tcp init method '%s'", init_method);
    return false;
  }

  host = endpoint.substr(0, colon);
  port = std::atoi(endpoint.c_str() + colon + 1);
  if (host.empty() || port <= 0) {
    set_error("invalid commux tcp init method '%s'", init_method);
    return false;
  }
  return true;
}

void set_default_env(const char* name, const char* value) {
  if (std::getenv(name) == nullptr) setenv(name, value, /*overwrite=*/0);
}

c10d::ReduceOp reduce_op(PNC_CommuxReduceOp op) {
  switch (op) {
    case PNC_COMMUX_SUM:
      return c10d::ReduceOp::SUM;
    case PNC_COMMUX_MIN:
      return c10d::ReduceOp::MIN;
    case PNC_COMMUX_MAX:
      return c10d::ReduceOp::MAX;
  }
  return c10d::ReduceOp::SUM;
}

Runtime* runtime_or_error() {
  if (!g_runtime || !g_runtime->pg) {
    set_error("commux runtime is not initialized");
    return nullptr;
  }
  return g_runtime.get();
}

torch::Tensor scalar_tensor(int64_t value) {
  return torch::full({1}, value,
                     torch::TensorOptions().dtype(torch::kInt64).device(
                         torch::kCPU));
}

int collective_error(Runtime* runtime, int local_status) {
  try {
    auto tensor = scalar_tensor(local_status == 0 ? 0 : 1);
    std::vector<at::Tensor> tensors{tensor};
    c10d::AllreduceOptions opts;
    opts.reduceOp = c10d::ReduceOp::SUM;
    runtime->pg->allreduce(tensors, opts)->wait();
    return tensor.item<int64_t>() == 0 ? PNC_COMMUX_SUCCESS
                                       : PNC_COMMUX_ERR_RUNTIME;
  } catch (const std::exception& e) {
    set_error("commux collective error check failed: %s", e.what());
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

}  // namespace

extern "C" {

const char* pnc_commux_last_error(void) { return g_last_error; }

int pnc_commux_init(const char* backend, const char* init_method, int rank,
                    int world_size) {
  std::lock_guard<std::mutex> lock(g_mutex);
  std::string host;
  int port = 0;

  if (backend != nullptr && backend[0] != '\0' &&
      std::strcmp(backend, "ucx") != 0 && std::strcmp(backend, "commux") != 0) {
    set_error("unsupported commux backend '%s'", backend);
    return PNC_COMMUX_ERR_ARG;
  }

  try {
    if (rank < 0) {
      const char* env_rank = getenv_required("RANK");
      if (env_rank == nullptr) return PNC_COMMUX_ERR_ARG;
      rank = std::atoi(env_rank);
    }
    if (world_size < 0) {
      const char* env_world = getenv_required("WORLD_SIZE");
      if (env_world == nullptr) return PNC_COMMUX_ERR_ARG;
      world_size = std::atoi(env_world);
    }
    if (rank < 0 || world_size <= 0 || rank >= world_size) {
      set_error("invalid commux rank/world_size rank=%d world_size=%d", rank,
                world_size);
      return PNC_COMMUX_ERR_ARG;
    }
    if (!parse_tcp_init_method(init_method, host, port))
      return PNC_COMMUX_ERR_ARG;

    set_default_env("UCX_TLS", "^cuda_copy,cuda_ipc,gdr_copy");

    c10d::TCPStoreOptions store_opts;
    store_opts.port = static_cast<std::uint16_t>(port);
    store_opts.numWorkers = world_size;
    store_opts.isServer = rank == 0;
    store_opts.waitWorkers = true;

    auto runtime = std::make_unique<Runtime>();
    runtime->store = c10::make_intrusive<c10d::TCPStore>(host, store_opts);
    runtime->pg = c10::make_intrusive<commux::ProcessGroupUCX>(
        runtime->store, rank, world_size);
    runtime->rank = rank;
    runtime->size = world_size;
    runtime->pg->barrier()->wait();

    g_runtime = std::move(runtime);
    g_last_error[0] = '\0';
    return PNC_COMMUX_SUCCESS;
  } catch (const std::exception& e) {
    set_error("commux init failed: %s", e.what());
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

int pnc_commux_init_env(const char* backend) {
  return pnc_commux_init(backend, "env://", -1, -1);
}

int pnc_commux_finalize(void) {
  std::lock_guard<std::mutex> lock(g_mutex);
  try {
    if (g_runtime && g_runtime->pg)
      g_runtime->pg->barrier()->wait();
    g_runtime.reset();
    g_last_error[0] = '\0';
    return PNC_COMMUX_SUCCESS;
  } catch (const std::exception& e) {
    set_error("commux finalize failed: %s", e.what());
    g_runtime.reset();
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

int pnc_commux_rank(int* rankp) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (rankp == nullptr) {
    set_error("pnc_commux_rank requires a non-null output pointer");
    return PNC_COMMUX_ERR_ARG;
  }
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  *rankp = runtime->rank;
  return PNC_COMMUX_SUCCESS;
}

int pnc_commux_size(int* sizep) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (sizep == nullptr) {
    set_error("pnc_commux_size requires a non-null output pointer");
    return PNC_COMMUX_ERR_ARG;
  }
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  *sizep = runtime->size;
  return PNC_COMMUX_SUCCESS;
}

int pnc_commux_barrier(void) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  try {
    runtime->pg->barrier()->wait();
    return PNC_COMMUX_SUCCESS;
  } catch (const std::exception& e) {
    set_error("commux barrier failed: %s", e.what());
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

int pnc_commux_allreduce_int64(int64_t value, PNC_CommuxReduceOp op,
                               int64_t* resultp) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (resultp == nullptr) {
    set_error("pnc_commux_allreduce_int64 requires a non-null result");
    return PNC_COMMUX_ERR_ARG;
  }
  if (op < PNC_COMMUX_SUM || op > PNC_COMMUX_MAX) {
    set_error("invalid commux reduce op %d", static_cast<int>(op));
    return PNC_COMMUX_ERR_ARG;
  }
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;

  try {
    auto tensor = scalar_tensor(value);
    std::vector<at::Tensor> tensors{tensor};
    c10d::AllreduceOptions opts;
    opts.reduceOp = reduce_op(op);
    runtime->pg->allreduce(tensors, opts)->wait();
    *resultp = tensor.item<int64_t>();
    return PNC_COMMUX_SUCCESS;
  } catch (const std::exception& e) {
    set_error("commux allreduce_int64 failed: %s", e.what());
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

int pnc_commux_send_int64(int64_t value, int dst, int tag) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  try {
    auto tensor = scalar_tensor(value);
    std::vector<at::Tensor> tensors{tensor};
    runtime->pg->send(tensors, dst, tag)->wait();
    return PNC_COMMUX_SUCCESS;
  } catch (const std::exception& e) {
    set_error("commux send_int64 failed: %s", e.what());
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

int pnc_commux_recv_int64(int src, int tag, int64_t* valuep) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (valuep == nullptr) {
    set_error("pnc_commux_recv_int64 requires a non-null result");
    return PNC_COMMUX_ERR_ARG;
  }
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  try {
    auto tensor = torch::zeros({1}, torch::TensorOptions().dtype(torch::kInt64));
    std::vector<at::Tensor> tensors{tensor};
    runtime->pg->recv(tensors, src, tag)->wait();
    *valuep = tensor.item<int64_t>();
    return PNC_COMMUX_SUCCESS;
  } catch (const std::exception& e) {
    set_error("commux recv_int64 failed: %s", e.what());
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

int pnc_commux_send_bytes(const void* data, int64_t nbytes, int dst, int tag) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  if (nbytes < 0 || (nbytes > 0 && data == nullptr)) {
    set_error("pnc_commux_send_bytes received invalid buffer");
    return PNC_COMMUX_ERR_ARG;
  }
  try {
    auto tensor = torch::empty({nbytes}, torch::TensorOptions().dtype(torch::kUInt8));
    if (nbytes > 0) std::memcpy(tensor.data_ptr(), data, static_cast<size_t>(nbytes));
    std::vector<at::Tensor> tensors{tensor};
    runtime->pg->send(tensors, dst, tag)->wait();
    return PNC_COMMUX_SUCCESS;
  } catch (const std::exception& e) {
    set_error("commux send_bytes failed: %s", e.what());
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

int pnc_commux_recv_bytes(void* data, int64_t nbytes, int src, int tag) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  if (nbytes < 0 || (nbytes > 0 && data == nullptr)) {
    set_error("pnc_commux_recv_bytes received invalid buffer");
    return PNC_COMMUX_ERR_ARG;
  }
  try {
    auto tensor = torch::empty({nbytes}, torch::TensorOptions().dtype(torch::kUInt8));
    std::vector<at::Tensor> tensors{tensor};
    runtime->pg->recv(tensors, src, tag)->wait();
    if (nbytes > 0) std::memcpy(data, tensor.data_ptr(), static_cast<size_t>(nbytes));
    return PNC_COMMUX_SUCCESS;
  } catch (const std::exception& e) {
    set_error("commux recv_bytes failed: %s", e.what());
    return PNC_COMMUX_ERR_RUNTIME;
  }
}

int pnc_commux_file_open(const char* path, int flags, int mode,
                         int64_t* handlep) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (handlep == nullptr || path == nullptr) {
    set_error("pnc_commux_file_open requires path and handle output");
    return PNC_COMMUX_ERR_ARG;
  }
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;

  int fd = ::open(path, flags, static_cast<mode_t>(mode));
  int local_status = fd >= 0 ? 0 : errno;
  if (local_status != 0)
    set_error("open('%s') failed on rank %d: %s", path, runtime->rank,
              std::strerror(local_status));
  int cerr = collective_error(runtime, local_status);
  if (cerr != PNC_COMMUX_SUCCESS) {
    if (fd >= 0) ::close(fd);
    return cerr;
  }

  int64_t handle = runtime->next_file_handle++;
  runtime->files[handle] = fd;
  *handlep = handle;
  return PNC_COMMUX_SUCCESS;
}

int pnc_commux_file_close(int64_t handle) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  auto it = runtime->files.find(handle);
  if (it == runtime->files.end()) {
    set_error("invalid commux file handle %" PRId64, handle);
    return PNC_COMMUX_ERR_ARG;
  }
  int fd = it->second;
  runtime->files.erase(it);
  int local_status = ::close(fd) == 0 ? 0 : errno;
  if (local_status != 0)
    set_error("close failed on rank %d: %s", runtime->rank,
              std::strerror(local_status));
  return collective_error(runtime, local_status);
}

int pnc_commux_file_pwrite(int64_t handle, const void* data, int64_t nbytes,
                           int64_t offset) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  if (nbytes < 0 || offset < 0 || (nbytes > 0 && data == nullptr)) {
    set_error("invalid pwrite arguments");
    return PNC_COMMUX_ERR_ARG;
  }
  auto it = runtime->files.find(handle);
  if (it == runtime->files.end()) {
    set_error("invalid commux file handle %" PRId64, handle);
    return PNC_COMMUX_ERR_ARG;
  }

  const char* ptr = static_cast<const char*>(data);
  int64_t done = 0;
  while (done < nbytes) {
    ssize_t n = ::pwrite(it->second, ptr + done,
                         static_cast<size_t>(nbytes - done),
                         static_cast<off_t>(offset + done));
    if (n < 0) {
      set_error("pwrite failed on rank %d: %s", runtime->rank,
                std::strerror(errno));
      return PNC_COMMUX_ERR_RUNTIME;
    }
    if (n == 0) {
      set_error("pwrite made no progress on rank %d", runtime->rank);
      return PNC_COMMUX_ERR_RUNTIME;
    }
    done += n;
  }
  return PNC_COMMUX_SUCCESS;
}

int pnc_commux_file_sync(int64_t handle) {
  std::lock_guard<std::mutex> lock(g_mutex);
  Runtime* runtime = runtime_or_error();
  if (runtime == nullptr) return PNC_COMMUX_ERR_RUNTIME;
  auto it = runtime->files.find(handle);
  if (it == runtime->files.end()) {
    set_error("invalid commux file handle %" PRId64, handle);
    return PNC_COMMUX_ERR_ARG;
  }
  int local_status = ::fsync(it->second) == 0 ? 0 : errno;
  if (local_status != 0)
    set_error("fsync failed on rank %d: %s", runtime->rank,
              std::strerror(local_status));
  return collective_error(runtime, local_status);
}

}  // extern "C"
