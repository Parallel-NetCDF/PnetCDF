#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2815 2017-01-01 23:06:33Z wkliao $
#
# @configure_input@

M4_HFILES = ncx_h.m4

M4FLAGS += -I${top_srcdir}/m4

ncx.h: ncx_h.m4 Makefile
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

# automake says "... BUILT_SOURCES is honored only by 'make all', 'make check',
# and 'make install'. This means you cannot build a specific target (e.g.,
# 'make foo') in a clean tree if it depends on a built source."
BUILT_SOURCES = ncx.h

EXTRA_DIST = common.h $(M4_HFILES)

CLEANFILES = ncx.h

tests-local: all

