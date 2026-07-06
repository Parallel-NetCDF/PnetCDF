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
    int ncid = -1, dim_rank = -1, dim_col = -1, varid = -1;
    int rank = atoi(getenv("RANK") ? getenv("RANK") : "0");
    int world = atoi(getenv("WORLD_SIZE") ? getenv("WORLD_SIZE") : "1");
    PNC_Offset start[2], count[2];
    float values[2];
    int req = -1, status = NC_NOERR;

    nerrs += check(world == 4, "expected four ranks");
    remove("ncmpix-cdf5-4proc.nc");

    nerrs += check_call(ncmpix_create(PNC_COMM_WORLD, "ncmpix-cdf5-4proc.nc",
                                      NC_CLOBBER, PNC_INFO_NULL, &ncid),
                        "ncmpix_create");
    nerrs += check_call(ncmpix_def_dim(ncid, "rank", 4, &dim_rank),
                        "ncmpix_def_dim(rank)");
    nerrs += check_call(ncmpix_def_dim(ncid, "col", 2, &dim_col),
                        "ncmpix_def_dim(col)");
    {
        int dims[2] = {dim_rank, dim_col};
        nerrs += check_call(ncmpix_def_var(ncid, "value", NC_FLOAT, 2, dims,
                                           &varid),
                            "ncmpix_def_var");
    }
    nerrs += check_call(ncmpix_put_att_text(ncid, varid, "units", 1, "1"),
                        "ncmpix_put_att_text");
    nerrs += check_call(ncmpix_enddef(ncid), "ncmpix_enddef");

    start[0] = rank;
    start[1] = 0;
    count[0] = 1;
    count[1] = 2;
    values[0] = (float)(10 + rank);
    values[1] = (float)(20 + rank);
    nerrs += check_call(ncmpix_iput_vara_float(ncid, varid, start, count,
                                               values, &req),
                        "ncmpix_iput_vara_float");
    nerrs += check_call(ncmpix_wait_all(ncid, 1, &req, &status),
                        "ncmpix_wait_all");
    nerrs += check(status == NC_NOERR, "request status should be NC_NOERR");
    nerrs += check_call(ncmpix_close(ncid), "ncmpix_close");

    if (rank == 0) {
        int rnc = -1, rvar = -1;
        float got[8];
        nerrs += check_call(nc_open("ncmpix-cdf5-4proc.nc", NC_NOWRITE, &rnc),
                            "nc_open verify");
        nerrs += check_call(nc_inq_varid(rnc, "value", &rvar),
                            "nc_inq_varid verify");
        nerrs += check_call(nc_get_var_float(rnc, rvar, got),
                            "nc_get_var_float verify");
        for (int r = 0; r < 4; ++r) {
            nerrs += check(got[2 * r] == (float)(10 + r),
                           "first column value mismatch");
            nerrs += check(got[2 * r + 1] == (float)(20 + r),
                           "second column value mismatch");
        }
        nc_close(rnc);
    }

    if (nerrs == 0 && rank == 0)
        printf("PASS: ncmpix CDF-5 4-process test\n");
    return nerrs == 0 ? 0 : 1;
}
