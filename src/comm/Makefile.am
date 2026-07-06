#
# Copyright (C) 2026, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

AM_CPPFLAGS  = -I${top_srcdir}/src/include
AM_CPPFLAGS += -I${top_builddir}/src/include
AM_CPPFLAGS += @COMMUX_CPPFLAGS@

AM_CFLAGS = @C_VISIBILITY_FLAGS@

noinst_LTLIBRARIES = libpnccomm.la

if COMM_BACKEND_COMMUX
libpnccomm_la_SOURCES = pnc_comm_commux.c
else
libpnccomm_la_SOURCES = pnc_comm_mpi.c
endif

EXTRA_DIST = pnc_comm_mpi.c pnc_comm_commux.c pnc_commux_runtime.cpp

tests-local: all
