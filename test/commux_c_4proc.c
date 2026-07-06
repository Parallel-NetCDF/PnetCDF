/*
 * C caller test for the Python-backed commux C ABI.
 */

#include <pnetcdf_commux.h>

#include <stdint.h>
#include <stdio.h>

static int
check(int ok, const char *msg)
{
    if (!ok) {
        fprintf(stderr, "FAIL: %s\n", msg);
        return 1;
    }
    return 0;
}

static int
check_call(int err, const char *call)
{
    if (err != PNC_COMMUX_SUCCESS) {
        fprintf(stderr, "FAIL: %s: %s\n", call, pnc_commux_last_error());
        return 1;
    }
    return 0;
}

int
main(void)
{
    int nerrs = 0;
    int rank = -1, size = -1;
    int64_t sum = 0;

    nerrs += check_call(pnc_commux_init_env("ucx"), "pnc_commux_init_env");
    nerrs += check_call(pnc_commux_rank(&rank), "pnc_commux_rank");
    nerrs += check_call(pnc_commux_size(&size), "pnc_commux_size");
    nerrs += check(size == 4, "expected four commux ranks");

    nerrs += check_call(pnc_commux_allreduce_int64((int64_t)rank + 1,
                                                   PNC_COMMUX_SUM, &sum),
                        "pnc_commux_allreduce_int64");
    nerrs += check(sum == 10, "expected allreduce sum of ranks 1..4");

    if (rank == 0) {
        nerrs += check_call(pnc_commux_send_int64(123, 1, 77),
                            "pnc_commux_send_int64");
    } else if (rank == 1) {
        int64_t value = 0;
        nerrs += check_call(pnc_commux_recv_int64(0, 77, &value),
                            "pnc_commux_recv_int64");
        nerrs += check(value == 123, "expected tagged value from rank 0");
    }

    nerrs += check_call(pnc_commux_barrier(), "pnc_commux_barrier");
    printf("rank %d/4 C commux wrapper ok\n", rank);
    nerrs += check_call(pnc_commux_finalize(), "pnc_commux_finalize");

    return nerrs == 0 ? 0 : 1;
}
