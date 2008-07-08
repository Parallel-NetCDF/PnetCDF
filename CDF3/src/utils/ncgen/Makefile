# Makefile for ncmpigen(1).
#
# $Id: Makefile.in 418 2004-10-29 21:09:56Z robl $

srcdir = .

INSTALL = /usr/bin/install -c

include ../../../macros.make

NCDUMP		= ../ncdump/ncmpidump
PROGRAM		= ncmpigen
INCLUDES	= -I$(srcdir)/../../../src/lib -I$(srcdir) -I.
MANUAL		= ncmpigen.1
GARBAGE		= $(PROGRAM) \
		  pnetcdf.inc \
		  c0.nc c1.cdl c1.nc c2.cdl \
		  f0.nc \
		  ctest.c ctest ctest0.nc ctest1.cdl \
		  ftest.f ftest ftest0.nc ftest1.cdl

lib_netcdf	= ../../../src/lib/libpnetcdf.a
ld_netcdf	= -L../../../src/lib -lpnetcdf

PACKING_LIST	= Makefile.in depend escapes.c generic.h \
		  ncgentab.c ncgentab.h \
		  genlib.c genlib.h getfill.c init.c load.c \
		  main.c ncmpigen.1 ncgen.h ncgen.l ncgen.y c0.cdl

OBJS		= main.o load.o ncgentab.o escapes.o \
		  getfill.o init.o genlib.o

lex		= flex
yacc		= bison -y

program_srcs	= main.c load.c ncgentab.c escapes.c \
		  getfill.c init.c genlib.c

all:		$(PROGRAM)

#  generating the fortran does not work yet
#test:           $(PROGRAM) b-test c-test f-test FORCE
test:           $(PROGRAM) b-test c-test FORCE

install:	
	$(INSTALL) -d -m 755 $(MANDIR)/man1
	$(INSTALL) -m 644 $(srcdir)/$(MANUAL) $(MANDIR)/man1/$(MANUAL)

	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 755 $(PROGRAM) $(BINDIR)/$(PROGRAM)

uninstall:
	-rm -f $(BINDIR)/$(PROGRAM)
	-rm -f $(MANDIR)/man1/$(MANUAL)

$(PROGRAM):	$(OBJS) $(lib_netcdf)
	$(LINK.c) $(OBJS) $(ld_netcdf) $(LIBS) 

ncgentab.c \
ncgentab.h:	ncgen.y ncgenyy.c ncgen.h
	$(yacc) -d $(srcdir)/ncgen.y; \
	mv y.tab.c ncgentab.c; \
	mv y.tab.h ncgentab.h

ncgenyy.c:	ncgen.l
	$(lex) $(srcdir)/ncgen.l; \
	mv lex.yy.c ncgenyy.c

vmstab.h \
vmstab.c:	ncgen.y
	@echo 1>&2 "$@ is out-of-date with respect to $?"
	@echo 1>&2 "It must be recreated via POSIX yacc(1) on a VMS system"
	false
vms_yy.c:	ncgenyy.c
	@echo 1>&2 "$@ is out-of-date with respect to $?"
	@echo 1>&2 "It must be recreated via POSIX lex(1) on a VMS system"
	false

#
# test "-b" option of ncgen
#
b-test:		$(PROGRAM) c1.cdl
	@./$(PROGRAM) -b c1.cdl && \
	$(NCDUMP) `pwd`/c1.nc > c2.cdl
	@if diff c1.cdl c2.cdl; then \
	    echo "*** $(PROGRAM) -b test successful ***"; \
	else \
	    echo "*** $(PROGRAM) -b test failed ***"; \
	    exit 1; \
	fi

#
# test "-c" option of ncgen
#
c-test:	$(PROGRAM) c1.cdl
	./$(PROGRAM) -c -o ctest0.nc $(srcdir)/c0.cdl > ctest.c && \
	$(CC) -o ctest $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) ctest.c $(ld_netcdf) $(LIBS) && \
	./ctest	&& \
	$(NCDUMP) -n c1 `pwd`/ctest0.nc > ctest1.cdl
	@if diff c1.cdl ctest1.cdl; then \
	    echo "*** $(PROGRAM) -c test successful ***"; \
	else \
	    echo "*** $(PROGRAM) -c test failed ***"; \
	    exit 1; \
	fi

c1.cdl:	$(PROGRAM) c0.cdl
	./$(PROGRAM) -b -o c0.nc $(srcdir)/c0.cdl
	$(NCDUMP) -n c1 `pwd`/c0.nc > $@

#
# test "-f" option of ncgen
#
f-test:	$(PROGRAM) c0.cdl c1.cdl
	@if [ -n "$(FC)" ]; then \
	    $(MAKE) $(MFLAGS) ftest1.cdl && \
	    if diff c1.cdl ftest1.cdl; then \
		echo "*** $(PROGRAM) -f test successful ***"; \
	    else \
		echo "*** $(PROGRAM) -f test failed ***"; \
		exit 1; \
	    fi; \
	else \
	    echo 1>&2 "\`$@' not made because no FORTRAN compiler"; \
	fi

ftest1.cdl: $(PROGRAM) c0.cdl pnetcdf.inc
	./$(PROGRAM) -f -o ftest0.nc $(srcdir)/c0.cdl > ftest.f
	$(COMPILE.f) ftest.f
	$(FC) $(FFLAGS) $(LDFLAGS) -o ftest ftest.o $(ld_netcdf) $(LIBS)
	./ftest
	$(NCDUMP) -n c1 ftest0.nc > ftest1.cdl

pnetcdf.inc:
	@if [ -n "$(FC)" ]; then \
	    cp $(srcdir)/../../libf/$@ .; \
	else \
	    echo 1>&2 "\`$@' not made because no FORTRAN compiler"; \
	fi

$(PROGRAM)_src : $(program_srcs)
	#setopt primary_language C
	#load -C $(CPPFLAGS) $(program_srcs)
	#load -C $(LIBS)
	#load -C /usr/lang/SC2.0.1/libansi.a
	#setopt program_name gribtonc

$(PROGRAM)_obj : $(program_srcs)
	#setopt primary_language C
	#load -C $(CPPFLAGS) $(OBJS)
	#load -C $(LIBS)
	#setopt program_name gribtonc

include $(srcdir)/../../../rules.make

include $(srcdir)/depend
