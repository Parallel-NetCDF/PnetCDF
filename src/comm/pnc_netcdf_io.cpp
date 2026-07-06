/*
 * Minimal UCX-coordinated CDF-5 file API for snapy pnetcdf output.
 */

#include <pnetcdf_comm.h>
#include <pnetcdf_commux.h>

#include <netcdf.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <string>
#include <vector>

namespace {

struct FileState {
  int ncid = -1;
  int rank = 0;
  int size = 1;
  int64_t file_handle = -1;
  int next_dimid = 0;
  int next_varid = 0;
  int next_reqid = 1;
  std::map<int, int> req_status;
  std::string path;
  bool define_mode = false;
  bool direct_io = false;
  PNC_Offset numrecs = 0;
  PNC_Offset record_size = 0;
  int unlimited_dimid = -1;
  struct VarLayout {
    int xtype = 0;
    PNC_Offset vsize = 0;
    PNC_Offset begin = 0;
    bool record = false;
    std::vector<PNC_Offset> dim_lengths;
  };
  std::map<int, VarLayout> vars;
};

std::map<int, FileState> g_files;
int g_next_handle = 1000;

int ensure_comm(void) {
  int rank = 0;
  int err = pnc_commux_rank(&rank);
  if (err == PNC_COMMUX_SUCCESS) return NC_NOERR;
  if (std::getenv("RANK") == nullptr && std::getenv("WORLD_SIZE") == nullptr)
    return NC_NOERR;
  err = pnc_commux_init_env("ucx");
  if (err == PNC_COMMUX_SUCCESS) return NC_NOERR;
  return NC_EINVAL;
}

int rank(void) {
  int r = 0;
  if (pnc_commux_rank(&r) != PNC_COMMUX_SUCCESS) return 0;
  return r;
}

int size(void) {
  int s = 1;
  if (pnc_commux_size(&s) != PNC_COMMUX_SUCCESS) return 1;
  return s;
}

int barrier(void) {
  if (size() <= 1) return NC_NOERR;
  return pnc_commux_barrier() == PNC_COMMUX_SUCCESS ? NC_NOERR : NC_EINVAL;
}

FileState* file_state(int ncid) {
  auto it = g_files.find(ncid);
  if (it == g_files.end()) return nullptr;
  return &it->second;
}

std::vector<size_t> to_size_t(const PNC_Offset* values, int n) {
  std::vector<size_t> out(n);
  for (int i = 0; i < n; ++i) out[i] = static_cast<size_t>(values[i]);
  return out;
}

int var_ndims(int ncid, int varid) {
  int ndims = 0;
  if (nc_inq_varndims(ncid, varid, &ndims) != NC_NOERR) return 0;
  return ndims;
}

PNC_Offset pad4(PNC_Offset n) { return (n + 3) & ~static_cast<PNC_Offset>(3); }

int xtype_size(int xtype) {
  switch (xtype) {
    case NC_BYTE:
    case NC_CHAR:
    case NC_UBYTE:
      return 1;
    case NC_SHORT:
    case NC_USHORT:
      return 2;
    case NC_INT:
    case NC_UINT:
    case NC_FLOAT:
      return 4;
    case NC_DOUBLE:
    case NC_INT64:
    case NC_UINT64:
      return 8;
    default:
      return 0;
  }
}

struct Cdf5Reader {
  std::vector<unsigned char> bytes;
  size_t pos = 0;
  bool ok = true;

  uint32_t u32() {
    if (pos + 4 > bytes.size()) {
      ok = false;
      return 0;
    }
    uint32_t v = (static_cast<uint32_t>(bytes[pos]) << 24) |
                 (static_cast<uint32_t>(bytes[pos + 1]) << 16) |
                 (static_cast<uint32_t>(bytes[pos + 2]) << 8) |
                 static_cast<uint32_t>(bytes[pos + 3]);
    pos += 4;
    return v;
  }

  uint64_t u64() {
    uint64_t v = 0;
    for (int i = 0; i < 8; ++i) v = (v << 8) | u32_byte();
    return v;
  }

  unsigned char u32_byte() {
    if (pos >= bytes.size()) {
      ok = false;
      return 0;
    }
    return bytes[pos++];
  }

  std::string name() {
    uint64_t len = u64();
    if (!ok || len > bytes.size() || pos + len > bytes.size()) {
      ok = false;
      return {};
    }
    std::string out(reinterpret_cast<const char*>(bytes.data() + pos),
                    static_cast<size_t>(len));
    pos += static_cast<size_t>(pad4(static_cast<PNC_Offset>(len)));
    if (pos > bytes.size()) ok = false;
    return out;
  }

  void skip(PNC_Offset n) {
    if (n < 0 || pos + static_cast<size_t>(n) > bytes.size()) {
      ok = false;
      return;
    }
    pos += static_cast<size_t>(n);
  }
};

bool parse_cdf5_layout(FileState& st) {
  std::ifstream in(st.path, std::ios::binary);
  if (!in) return false;
  Cdf5Reader r;
  r.bytes.assign(std::istreambuf_iterator<char>(in),
                 std::istreambuf_iterator<char>());
  if (r.bytes.size() < 12 || r.bytes[0] != 'C' || r.bytes[1] != 'D' ||
      r.bytes[2] != 'F' || r.bytes[3] != 5)
    return false;

  r.pos = 4;
  st.numrecs = static_cast<PNC_Offset>(r.u64());

  uint32_t tag = r.u32();
  uint64_t ndims = r.u64();
  if (!r.ok || (tag != 0 && tag != 10)) return false;
  std::vector<PNC_Offset> dim_lengths;
  for (uint64_t i = 0; i < ndims; ++i) {
    (void)r.name();
    PNC_Offset len = static_cast<PNC_Offset>(r.u64());
    if (!r.ok) return false;
    if (len == 0 && st.unlimited_dimid < 0) st.unlimited_dimid = i;
    dim_lengths.push_back(len);
  }

  tag = r.u32();
  uint64_t ngatts = r.u64();
  if (!r.ok || (tag != 0 && tag != 12)) return false;
  for (uint64_t i = 0; i < ngatts; ++i) {
    (void)r.name();
    int xtype = static_cast<int>(r.u32());
    uint64_t nelems = r.u64();
    int size = xtype_size(xtype);
    if (size <= 0) return false;
    r.skip(pad4(static_cast<PNC_Offset>(size) *
                static_cast<PNC_Offset>(nelems)));
  }

  tag = r.u32();
  uint64_t nvars = r.u64();
  if (!r.ok || (tag != 0 && tag != 11)) return false;
  st.vars.clear();
  st.record_size = 0;
  for (uint64_t varid = 0; varid < nvars; ++varid) {
    (void)r.name();
    uint64_t var_ndims = r.u64();
    std::vector<uint64_t> dimids;
    for (uint64_t i = 0; i < var_ndims; ++i) dimids.push_back(r.u64());

    tag = r.u32();
    uint64_t natts = r.u64();
    if (!r.ok || (tag != 0 && tag != 12)) return false;
    for (uint64_t i = 0; i < natts; ++i) {
      (void)r.name();
      int xtype = static_cast<int>(r.u32());
      uint64_t nelems = r.u64();
      int size = xtype_size(xtype);
      if (size <= 0) return false;
      r.skip(pad4(static_cast<PNC_Offset>(size) *
                  static_cast<PNC_Offset>(nelems)));
    }

    FileState::VarLayout layout;
    layout.xtype = static_cast<int>(r.u32());
    layout.vsize = static_cast<PNC_Offset>(r.u64());
    layout.begin = static_cast<PNC_Offset>(r.u64());
    if (!r.ok) return false;
    layout.record =
        !dimids.empty() && static_cast<int>(dimids[0]) == st.unlimited_dimid;
    for (uint64_t dimid : dimids) {
      if (dimid >= dim_lengths.size()) return false;
      layout.dim_lengths.push_back(dim_lengths[static_cast<size_t>(dimid)]);
    }
    st.vars[static_cast<int>(varid)] = layout;
    if (layout.record) st.record_size += layout.vsize;
  }
  return r.ok;
}

void store_be64(unsigned char out[8], PNC_Offset value) {
  uint64_t v = static_cast<uint64_t>(value);
  for (int i = 7; i >= 0; --i) {
    out[i] = static_cast<unsigned char>(v & 0xff);
    v >>= 8;
  }
}

void append_float_be(std::vector<unsigned char>& out, float value) {
  uint32_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  out.push_back(static_cast<unsigned char>((bits >> 24) & 0xff));
  out.push_back(static_cast<unsigned char>((bits >> 16) & 0xff));
  out.push_back(static_cast<unsigned char>((bits >> 8) & 0xff));
  out.push_back(static_cast<unsigned char>(bits & 0xff));
}

PNC_Offset row_major_offset(const std::vector<PNC_Offset>& dims,
                            const std::vector<PNC_Offset>& index) {
  PNC_Offset offset = 0;
  for (size_t i = 0; i < dims.size(); ++i)
    offset = offset * dims[i] + index[i];
  return offset;
}

int direct_write_float(FileState& st, int varid, const PNC_Offset* start,
                       const PNC_Offset* count, const float* value,
                       PNC_Offset* endrecp) {
  if (!st.direct_io || st.file_handle < 0) return NC_ENOTBUILT;
  auto it = st.vars.find(varid);
  if (it == st.vars.end()) return NC_ENOTBUILT;
  const auto& layout = it->second;
  if (layout.xtype != NC_FLOAT || layout.dim_lengths.empty())
    return NC_ENOTBUILT;

  int ndims = static_cast<int>(layout.dim_lengths.size());
  if (endrecp != nullptr && layout.record)
    *endrecp = std::max(*endrecp, start[0] + count[0]);

  int elem_size = xtype_size(layout.xtype);
  int last = ndims - 1;
  std::vector<PNC_Offset> prefix(std::max(0, last), 0);
  PNC_Offset value_base = 0;

  std::function<int(int, PNC_Offset)> write_runs =
      [&](int dim, PNC_Offset linear_prefix) -> int {
    if (dim == last) {
      std::vector<PNC_Offset> file_index(ndims);
      for (int i = 0; i < last; ++i) file_index[i] = start[i] + prefix[i];
      file_index[last] = start[last];

      PNC_Offset file_offset = 0;
      if (layout.record) {
        std::vector<PNC_Offset> rec_dims(layout.dim_lengths.begin() + 1,
                                         layout.dim_lengths.end());
        std::vector<PNC_Offset> rec_index(file_index.begin() + 1,
                                          file_index.end());
        file_offset = layout.begin + file_index[0] * st.record_size +
                      row_major_offset(rec_dims, rec_index) * elem_size;
      } else {
        file_offset =
            layout.begin +
            row_major_offset(layout.dim_lengths, file_index) * elem_size;
      }

      std::vector<unsigned char> bytes;
      bytes.reserve(static_cast<size_t>(count[last]) * sizeof(float));
      PNC_Offset value_offset = linear_prefix * count[last];
      for (PNC_Offset i = 0; i < count[last]; ++i)
        append_float_be(bytes, value[value_offset + i]);
      return pnc_commux_file_pwrite(
                 st.file_handle, bytes.data(),
                 static_cast<int64_t>(bytes.size()),
                 static_cast<int64_t>(file_offset)) == PNC_COMMUX_SUCCESS
                 ? NC_NOERR
                 : NC_EINVAL;
    }

    for (PNC_Offset i = 0; i < count[dim]; ++i) {
      prefix[dim] = i;
      int err = write_runs(dim + 1, linear_prefix * count[dim] + i);
      if (err != NC_NOERR) return err;
    }
    return NC_NOERR;
  };

  if (ndims == 0) return NC_ENOTBUILT;
  return write_runs(0, value_base);
}

int finish_direct_write(FileState& st, int local_status, PNC_Offset endrec) {
  int64_t failed = local_status == NC_NOERR ? 0 : 1;
  int64_t any_failed = 0;
  if (pnc_commux_allreduce_int64(failed, PNC_COMMUX_SUM, &any_failed) !=
      PNC_COMMUX_SUCCESS)
    return NC_EINVAL;
  int status = any_failed == 0 ? NC_NOERR : NC_EINVAL;

  int64_t local_endrec = static_cast<int64_t>(endrec);
  int64_t max_endrec = 0;
  if (pnc_commux_allreduce_int64(local_endrec, PNC_COMMUX_MAX, &max_endrec) !=
      PNC_COMMUX_SUCCESS)
    status = NC_EINVAL;
  if (status == NC_NOERR && max_endrec > st.numrecs) {
    unsigned char nrec[8];
    store_be64(nrec, static_cast<PNC_Offset>(max_endrec));
    if (pnc_commux_file_pwrite(st.file_handle, nrec, sizeof(nrec), 4) !=
        PNC_COMMUX_SUCCESS)
      status = NC_EINVAL;
    st.numrecs = static_cast<PNC_Offset>(max_endrec);
  }
  int berr = barrier();
  return status != NC_NOERR ? status : berr;
}

}  // namespace

extern "C" {

const char* ncmpix_strerror(int err) {
  if (err == NC_ENOTBUILT) return "feature was not built";
  return nc_strerror(err);
}

int ncmpix_create(PNC_Comm, const char* path, int cmode, PNC_Info, int* ncidp) {
  if (ncidp == nullptr) return NC_EINVAL;
  int err = ensure_comm();
  if (err != NC_NOERR) return err;

  FileState st;
  st.rank = rank();
  st.size = size();
  st.path = path != nullptr ? path : "";
  st.define_mode = true;
  if (st.path.empty()) return NC_EINVAL;

  if (st.rank == 0) {
    int mode = NC_CDF5;
    if (cmode & NC_NOCLOBBER)
      mode |= NC_NOCLOBBER;
    else
      mode |= NC_CLOBBER;
    err = nc_create(st.path.c_str(), mode, &st.ncid);
    if (err != NC_NOERR) return err;
    int old_fill = 0;
    nc_set_fill(st.ncid, NC_NOFILL, &old_fill);
  }
  if ((err = barrier()) != NC_NOERR) return err;

  int handle = g_next_handle++;
  g_files[handle] = st;
  *ncidp = handle;
  return NC_NOERR;
}

int ncmpix_open(PNC_Comm, const char* path, int omode, PNC_Info, int* ncidp) {
  if (ncidp == nullptr || path == nullptr) return NC_EINVAL;
  int err = ensure_comm();
  if (err != NC_NOERR) return err;

  FileState st;
  st.rank = rank();
  st.size = size();
  st.path = path;
  err = nc_open(path, omode, &st.ncid);
  if (err != NC_NOERR) return err;

  int handle = g_next_handle++;
  g_files[handle] = st;
  *ncidp = handle;
  return NC_NOERR;
}

int ncmpix_close(int ncid) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  int err = NC_NOERR;
  if (st->file_handle >= 0) {
    int serr = pnc_commux_file_sync(st->file_handle) == PNC_COMMUX_SUCCESS
                   ? NC_NOERR
                   : NC_EINVAL;
    int cerr = pnc_commux_file_close(st->file_handle) == PNC_COMMUX_SUCCESS
                   ? NC_NOERR
                   : NC_EINVAL;
    if (err == NC_NOERR) err = serr != NC_NOERR ? serr : cerr;
    st->file_handle = -1;
  }
  if (st->ncid >= 0) err = nc_close(st->ncid);
  g_files.erase(ncid);
  int berr = barrier();
  return err != NC_NOERR ? err : berr;
}

int ncmpix_def_dim(int ncid, const char* name, PNC_Offset len, int* idp) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  int id = st->next_dimid++;
  int err = NC_NOERR;
  if (st->rank == 0)
    err = nc_def_dim(st->ncid, name, static_cast<size_t>(len), &id);
  if (idp != nullptr) *idp = id;
  return err;
}

int ncmpix_def_var(int ncid, const char* name, nc_type xtype, int ndims,
                   const int* dimids, int* varidp) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  int id = st->next_varid++;
  int err = NC_NOERR;
  if (st->rank == 0)
    err = nc_def_var(st->ncid, name, xtype, ndims, dimids, &id);
  if (varidp != nullptr) *varidp = id;
  return err;
}

int ncmpix_put_att_text(int ncid, int varid, const char* name, PNC_Offset len,
                        const char* value) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  if (st->rank != 0) return NC_NOERR;
  return nc_put_att_text(st->ncid, varid, name, static_cast<size_t>(len), value);
}

int ncmpix_put_att_int(int ncid, int varid, const char* name, nc_type xtype,
                       PNC_Offset len, const int* value) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  if (st->rank != 0) return NC_NOERR;
  return nc_put_att_int(st->ncid, varid, name, xtype, static_cast<size_t>(len),
                        value);
}

int ncmpix_put_att_float(int ncid, int varid, const char* name, nc_type xtype,
                         PNC_Offset len, const float* value) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  if (st->rank != 0) return NC_NOERR;
  return nc_put_att_float(st->ncid, varid, name, xtype, static_cast<size_t>(len),
                          value);
}

int ncmpix_enddef(int ncid) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  int err = NC_NOERR;
  if (st->rank == 0) {
    err = nc_enddef(st->ncid);
    if (err == NC_NOERR) err = nc_close(st->ncid);
    st->ncid = -1;
  }
  if (err != NC_NOERR) return err;
  if ((err = barrier()) != NC_NOERR) return err;
  err = nc_open(st->path.c_str(), NC_WRITE, &st->ncid);
  if (err != NC_NOERR) return err;
  st->define_mode = false;
  st->direct_io = false;
  if (st->size > 1 && parse_cdf5_layout(*st)) {
    int64_t handle = -1;
    if (pnc_commux_file_open(st->path.c_str(), O_RDWR, 0666, &handle) ==
        PNC_COMMUX_SUCCESS) {
      st->file_handle = handle;
      st->direct_io = true;
    }
  }
  return barrier();
}

int ncmpix_put_vara_float_all(int ncid, int varid, const PNC_Offset* start,
                              const PNC_Offset* count, const float* value) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  PNC_Offset endrec = st->numrecs;
  int direct_status = direct_write_float(*st, varid, start, count, value,
                                         &endrec);
  if (direct_status != NC_ENOTBUILT)
    return finish_direct_write(*st, direct_status, endrec);

  int err = NC_NOERR;
  if (st->rank == 0) {
    int ndims = var_ndims(st->ncid, varid);
    auto s = to_size_t(start, ndims);
    auto c = to_size_t(count, ndims);
    err = nc_put_vara_float(st->ncid, varid, s.data(), c.data(), value);
  }
  int berr = barrier();
  return err != NC_NOERR ? err : berr;
}

int ncmpix_iput_vara_float(int ncid, int varid, const PNC_Offset* start,
                           const PNC_Offset* count, const float* value,
                           int* reqidp) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  int ndims = var_ndims(st->ncid, varid);
  auto s = to_size_t(start, ndims);
  auto c = to_size_t(count, ndims);
  int req = st->next_reqid++;
  int meta_tag = 9000 + req * 2;
  int data_tag = meta_tag + 1;
  PNC_Offset nvals = 1;
  for (int i = 0; i < ndims; ++i) nvals *= count[i];
  PNC_Offset nbytes = nvals * static_cast<PNC_Offset>(sizeof(float));
  int status = NC_NOERR;

  PNC_Offset endrec = st->numrecs;
  int direct_status = direct_write_float(*st, varid, start, count, value,
                                         &endrec);
  if (direct_status != NC_ENOTBUILT) {
    status = finish_direct_write(*st, direct_status, endrec);
    st->req_status[req] = status;
    if (reqidp != nullptr) *reqidp = req;
    return status;
  }

  if (st->rank == 0) {
    status = nc_put_vara_float(st->ncid, varid, s.data(), c.data(), value);
    for (int r = 1; r < st->size; ++r) {
      std::vector<PNC_Offset> meta(1 + 2 * ndims);
      if (pnc_commux_recv_bytes(meta.data(),
                                static_cast<int64_t>(meta.size() *
                                                     sizeof(PNC_Offset)),
                                r, meta_tag) != PNC_COMMUX_SUCCESS) {
        status = NC_EINVAL;
        continue;
      }
      int remote_ndims = static_cast<int>(meta[0]);
      if (remote_ndims != ndims) {
        status = NC_EINVAL;
        continue;
      }
      std::vector<size_t> rs(ndims), rc(ndims);
      PNC_Offset remote_vals = 1;
      for (int i = 0; i < ndims; ++i) {
        rs[i] = static_cast<size_t>(meta[1 + i]);
        rc[i] = static_cast<size_t>(meta[1 + ndims + i]);
        remote_vals *= meta[1 + ndims + i];
      }
      std::vector<float> remote(static_cast<size_t>(remote_vals));
      if (pnc_commux_recv_bytes(remote.data(),
                                static_cast<int64_t>(remote.size() *
                                                     sizeof(float)),
                                r, data_tag) != PNC_COMMUX_SUCCESS) {
        status = NC_EINVAL;
        continue;
      }
      int werr = nc_put_vara_float(st->ncid, varid, rs.data(), rc.data(),
                                   remote.data());
      if (status == NC_NOERR && werr != NC_NOERR) status = werr;
    }
  } else {
    std::vector<PNC_Offset> meta(1 + 2 * ndims);
    meta[0] = ndims;
    for (int i = 0; i < ndims; ++i) {
      meta[1 + i] = start[i];
      meta[1 + ndims + i] = count[i];
    }
    if (pnc_commux_send_bytes(meta.data(),
                              static_cast<int64_t>(meta.size() *
                                                   sizeof(PNC_Offset)),
                              0, meta_tag) != PNC_COMMUX_SUCCESS ||
        pnc_commux_send_bytes(value, static_cast<int64_t>(nbytes), 0,
                              data_tag) != PNC_COMMUX_SUCCESS) {
      status = NC_EINVAL;
    }
  }
  int berr = barrier();
  if (status == NC_NOERR && berr != NC_NOERR) status = berr;
  st->req_status[req] = status;
  if (reqidp != nullptr) *reqidp = req;
  return status;
}

int ncmpix_iput_var_float(int ncid, int varid, const float* value,
                          int* reqidp) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  int status = NC_NOERR;
  for (int r = 0; r < st->size; ++r) {
    if (st->rank == r) status = nc_put_var_float(st->ncid, varid, value);
    int berr = barrier();
    if (status == NC_NOERR && berr != NC_NOERR) status = berr;
  }
  int req = st->next_reqid++;
  st->req_status[req] = status;
  if (reqidp != nullptr) *reqidp = req;
  return status;
}

int ncmpix_wait_all(int ncid, int num, int* reqids, int* statuses) {
  FileState* st = file_state(ncid);
  if (st == nullptr) return NC_EBADID;
  int err = NC_NOERR;
  for (int i = 0; i < num; ++i) {
    int status = NC_NOERR;
    auto it = st->req_status.find(reqids[i]);
    if (it != st->req_status.end()) {
      status = it->second;
      st->req_status.erase(it);
    }
    if (statuses != nullptr) statuses[i] = status;
    if (err == NC_NOERR && status != NC_NOERR) err = status;
  }
  int berr = barrier();
  return err != NC_NOERR ? err : berr;
}

}  // extern "C"
