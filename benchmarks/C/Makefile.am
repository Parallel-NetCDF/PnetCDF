#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

srcdir  = @srcdir@
VPATH   = @srcdir@

include ../../macros.make

# note the order of -L list matters
INCLUDES  = -I../../src/lib
LDFLAGS  := -L../../src/lib $(LDFLAGS)
LIBS     := -lpnetcdf $(LIBS)

C_SRCS    = aggregation.c \
            write_block_read_column.c

PROGS     = $(C_SRCS:.c=)
OBJS      = $(C_SRCS:.c=.o)

GARBAGE      = $(PROGS) *.nc

PACKING_LIST = $(C_SRCS) Makefile.in depend

all: $(PROGS)

install:

uninstall:

aggregation: aggregation.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

write_block_read_column: write_block_read_column.o $(LIBRARY)
	$(LINK.c) $< $(LDFLAGS) $(LIBS)

TEST_MPIRUN_4  = $(subst NP,4,$(TEST_MPIRUN))

ptest: $(PROGS)
	@for i in $(PROGS); do { \
	$(TEST_MPIRUN_4) ./$$i -q 10 $(TEST_OUTDIR)/testfile.nc ; \
	if [ $$? = 0 ] ; then \
	    echo "PASS:  C  parallel run on 4 processes --------------- $$i"; \
	else \
	    echo "FAILED:  C  parallel run on 4 processes ------------- $$i"; \
	    exit 1; \
	fi ; } ; done

ptests:

include $(srcdir)/depend

include $(srcdir)/../../rules.make

$(LIBRARY): ;

