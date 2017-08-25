#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

srcdir = @srcdir@
VPATH  = @srcdir@

include ../../macros.make

INCLUDES  = -I../../src/lib -I$(srcdir)/../common
LDFLAGS  := $(LDFLAGS) -L../common
LIBS     := $(LIBRARY) -ltestutils $(LIBS)

SRCS    = test_read.c \
          test_read_indep.c \
          test_write.c \
          test_write_indep.c

OBJS    = $(SRCS:.c=.o)
PROGS   = $(SRCS:.c=)

GARBAGE      = $(PROGS)
PACKING_LIST = $(SRCS) Makefile.in README

all: $(PROGS)
	set -e; cd ../data && $(MAKE) $(MFLAGS) all

$(OBJS): $(srcdir)/../common/testutils.h

$(PROGS): ../common/libtestutils.a

../common/libtestutils.a:
	set -e; cd ../common && $(MAKE) $(MFLAGS) all

test_read: test_read.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

test_read_indep: test_read_indep.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

test_write: test_write.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

test_write_indep: test_write_indep.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

testing check verbose_testing:

ptest ptests ptest2 ptest4 ptest6 ptest8 ptest10:

include $(srcdir)/../../rules.make

$(LIBRARY): ;

