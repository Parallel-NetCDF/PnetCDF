#
# Copyright (C) 2015, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

bin_PROGRAMS = ncoffsets
EXTRA_DIST = ncoffsets.c

ncoffsets$(EXEEXT): ncoffsets.c
	$(SEQ_CC) -o $@ $<

dist_man_MANS = ncoffsets.1

tests-local: all

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

