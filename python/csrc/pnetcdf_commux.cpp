#include <pnetcdf_commux.h>

#include <pybind11/pybind11.h>

#include <stdexcept>

namespace py = pybind11;

namespace {

void check(int err) {
  if (err != PNC_COMMUX_SUCCESS)
    throw std::runtime_error(pnc_commux_last_error());
}

}  // namespace

PYBIND11_MODULE(_pnetcdf_commux, m) {
  m.doc() = "Thin Python binding for the PnetCDF commux C ABI";

  py::enum_<PNC_CommuxReduceOp>(m, "ReduceOp")
      .value("SUM", PNC_COMMUX_SUM)
      .value("MIN", PNC_COMMUX_MIN)
      .value("MAX", PNC_COMMUX_MAX);

  m.def("init", [](const char* backend, const char* init_method, int rank,
                   int world_size) {
    check(pnc_commux_init(backend, init_method, rank, world_size));
  });
  m.def("init_env", [](const char* backend) {
    check(pnc_commux_init_env(backend));
  }, py::arg("backend") = "ucx");
  m.def("finalize", []() { check(pnc_commux_finalize()); });
  m.def("rank", []() {
    int rank = -1;
    check(pnc_commux_rank(&rank));
    return rank;
  });
  m.def("size", []() {
    int size = -1;
    check(pnc_commux_size(&size));
    return size;
  });
  m.def("barrier", []() { check(pnc_commux_barrier()); });
  m.def("allreduce_int64", [](int64_t value, PNC_CommuxReduceOp op) {
    int64_t result = 0;
    check(pnc_commux_allreduce_int64(value, op, &result));
    return result;
  });
  m.def("send_int64", [](int64_t value, int dst, int tag) {
    check(pnc_commux_send_int64(value, dst, tag));
  });
  m.def("recv_int64", [](int src, int tag) {
    int64_t value = 0;
    check(pnc_commux_recv_int64(src, tag, &value));
    return value;
  });
}
