#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.am 2847 2017-01-17 07:40:16Z wkliao $
#
# @configure_input@

SUFFIXES = .o .c .h

noinst_LTLIBRARIES = libdispatchers.la

libdispatchers_la_SOURCES = file.c

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

