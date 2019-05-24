/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef TYPES_H
#define TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

enum BUFFERING_STATE { 
    buffering_stopped = 0,     /* Buffering stopped after overflow */
    buffering_ongoing = 1,     /* Buffering is going on (default value at start) */
};

enum BUFFERING_STRATEGY { 
    no_buffering            = 0, /* No buffering by common layer */
    stop_on_overflow        = 1, /* Buffering stopped after overflow. Does not demolish existing buffer until close()*/
    continue_with_new_pg    = 2  /* Buffering is going on with new PG starting after overflow */
    //contiguous_buffering    = 3 /* NOT SUPPORTED: Buffering is going on (default value at start) */
};

#ifdef __cplusplus
}
#endif

#endif
