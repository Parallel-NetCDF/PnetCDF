#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

srcdir  = @srcdir@
VPATH   = @srcdir@

include ../../macros.make

INCLUDES  = -I../../src/lib -I$(srcdir)/../common
LDFLAGS  := $(LDFLAGS) -L../common
LIBS     := $(LIBRARY) -ltestutils $(LIBS)

SRCS    = test_subarray.c \
          test_darray.c \
          test_nonblocking.c

HEADERS = test_dtype.h

OBJS    = $(SRCS:.c=.o)
PROGS   = $(SRCS:.c=)

GARBAGE      = $(PROGS)
PACKING_LIST = $(SRCS) $(HEADERS) Makefile.in README

all: $(PROGS)

$(OBJS): $(srcdir)/../common/testutils.h

$(PROGS): ../common/libtestutils.a

../common/libtestutils.a:
	set -e; cd ../common && $(MAKE) $(MFLAGS) all

test_subarray: test_subarray.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

test_darray: test_darray.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

test_nonblocking: test_nonblocking.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

testing check verbose_testing:

ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:

include $(srcdir)/../../rules.make

$(LIBRARY): ;

