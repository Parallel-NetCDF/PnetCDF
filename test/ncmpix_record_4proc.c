#include <pnetcdf_comm.h>

#include <netcdf.h>

#include <stdio.h>
#include <stdlib.h>

static int check(int ok, const char *msg) {
    if (!ok) {
        fprintf(stderr, "FAIL: %s\n", msg);
        return 1;
    }
    return 0;
}

static int check_call(int err, const char *call) {
    if (err != NC_NOERR) {
        fprintf(stderr, "FAIL: %s: %s\n", call, ncmpix_strerror(err));
        return 1;
    }
    return 0;
}

int main(void) {
    int nerrs = 0;
    int ncid = -1, dim_time = -1, dim_x = -1, time_id = -1, rec_id = -1;
    int rank = atoi(getenv("RANK") ? getenv("RANK") : "0");
    int world = atoi(getenv("WORLD_SIZE") ? getenv("WORLD_SIZE") : "1");
    int req = -1, status = NC_NOERR;
    PNC_Offset one_start[1] = {0}, one_count[1] = {1};
    PNC_Offset start[2], count[2];
    float time_value = 3.5f;
    float values[2];

    nerrs += check(world == 4, "expected four ranks");
    remove("ncmpix-record-4proc.nc");

    nerrs += check_call(ncmpix_create(PNC_COMM_WORLD, "ncmpix-record-4proc.nc",
                                      NC_CLOBBER, PNC_INFO_NULL, &ncid),
                        "ncmpix_create");
    nerrs += check_call(ncmpix_def_dim(ncid, "time", NC_UNLIMITED, &dim_time),
                        "ncmpix_def_dim(time)");
    nerrs += check_call(ncmpix_def_dim(ncid, "x", 8, &dim_x),
                        "ncmpix_def_dim(x)");
    nerrs += check_call(ncmpix_def_var(ncid, "time", NC_FLOAT, 1, &dim_time,
                                       &time_id),
                        "ncmpix_def_var(time)");
    {
        int dims[2] = {dim_time, dim_x};
        nerrs += check_call(ncmpix_def_var(ncid, "rec", NC_FLOAT, 2, dims,
                                           &rec_id),
                            "ncmpix_def_var(rec)");
    }
    nerrs += check_call(ncmpix_enddef(ncid), "ncmpix_enddef");

    nerrs += check_call(ncmpix_put_vara_float_all(ncid, time_id, one_start,
                                                  one_count, &time_value),
                        "ncmpix_put_vara_float_all(time)");

    start[0] = 0;
    start[1] = rank * 2;
    count[0] = 1;
    count[1] = 2;
    values[0] = (float)(100 + rank * 10);
    values[1] = (float)(101 + rank * 10);
    nerrs += check_call(ncmpix_iput_vara_float(ncid, rec_id, start, count,
                                               values, &req),
                        "ncmpix_iput_vara_float(rec)");
    nerrs += check_call(ncmpix_wait_all(ncid, 1, &req, &status),
                        "ncmpix_wait_all");
    nerrs += check(status == NC_NOERR, "request status should be NC_NOERR");
    nerrs += check_call(ncmpix_close(ncid), "ncmpix_close");

    if (rank == 0) {
        int rnc = -1, rtime = -1, rrec = -1;
        size_t ntime = 0;
        float got_time = 0.0f;
        float got[8];
        nerrs += check_call(nc_open("ncmpix-record-4proc.nc", NC_NOWRITE,
                                    &rnc),
                            "nc_open verify");
        nerrs += check_call(nc_inq_dimlen(rnc, dim_time, &ntime),
                            "nc_inq_dimlen(time)");
        nerrs += check(ntime == 1, "record length should be 1");
        nerrs += check_call(nc_inq_varid(rnc, "time", &rtime),
                            "nc_inq_varid(time)");
        nerrs += check_call(nc_inq_varid(rnc, "rec", &rrec),
                            "nc_inq_varid(rec)");
        nerrs += check_call(nc_get_var_float(rnc, rtime, &got_time),
                            "nc_get_var_float(time)");
        nerrs += check(got_time == time_value, "time value mismatch");
        nerrs += check_call(nc_get_var_float(rnc, rrec, got),
                            "nc_get_var_float(rec)");
        for (int r = 0; r < 4; ++r) {
            nerrs += check(got[2 * r] == (float)(100 + r * 10),
                           "record first value mismatch");
            nerrs += check(got[2 * r + 1] == (float)(101 + r * 10),
                           "record second value mismatch");
        }
        nc_close(rnc);
    }

    if (nerrs == 0 && rank == 0)
        printf("PASS: ncmpix record CDF-5 4-process test\n");
    return nerrs == 0 ? 0 : 1;
}
