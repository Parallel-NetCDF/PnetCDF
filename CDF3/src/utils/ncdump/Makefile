# Makefile for ncmpidump(1).
#
# $Id: Makefile.in 322 2003-12-12 17:38:08Z robl $

srcdir = .

INSTALL = /usr/bin/install -c

include ../../../macros.make

NCGEN		= ../ncgen/ncmpigen
INCLUDES	= -I$(srcdir)/../../../src/lib -I$(srcdir)

c_sources	= ncdump.c vardata.c dumplib.c
headers		= ncdump.h vardata.h dumplib.h

PROGRAM		= ncmpidump

PACKING_LIST	= $(c_sources) $(headers) depend \
		  test0.cdl ncmpidump.1 Makefile
MANUAL		= ncmpidump.1

lib_netcdf	= ../../../src/lib/libpnetcdf.a
ld_netcdf	= -L../../../src/lib -lpnetcdf

OBJS		=  ncdump.o vardata.o dumplib.o
GARBAGE		= $(PROGRAM) test0.nc test1.nc test1.cdl test2.cdl

all:		$(PROGRAM)

$(PROGRAM):	$(lib_netcdf) $(OBJS)
	$(LINK.c) $(OBJS) $(ld_netcdf) $(LIBS) 

test:		$(PROGRAM) FORCE
	$(NCGEN) -b $(srcdir)/test0.cdl -o `pwd`/test0.nc
	./$(PROGRAM) -n test1 test0.nc > test1.cdl
	$(NCGEN) -b `pwd`/test1.cdl -o `pwd`/test1.nc
	./$(PROGRAM) test1.nc > test2.cdl
	@cmp test1.cdl test2.cdl && \
	     echo "*** $(PROGRAM) test successful ***"

install:	
	$(INSTALL) -d -m 755 $(MANDIR)/manl
	$(INSTALL) -m 644 $(srcdir)/$(MANUAL) $(MANDIR)/man1/$(MANUAL)

	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 755 $(PROGRAM) $(BINDIR)/$(PROGRAM)

uninstall:
	-rm -f $(BINDIR)/$(PROGRAM)
	-rm -f $(MANDIR)/man1/$(MANUAL)

$(PROGRAM)_oc : $(c_sources)
	#setopt primary_language C
	#load -C $(CPPFLAGS) $(c_sources)
	#load -C $(LIBS)
	#setopt program_name $(PROGRAM)

TAGS:		FORCE
	etags `echo $(PACKING_LIST) | fmt -1 | egrep '\.c|\.h'

include $(srcdir)/../../../rules.make
include $(srcdir)/depend
