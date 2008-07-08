srcdir = .

SUBDIRS = ncgen ncdump

all: ncmpigen ncmpidump

ncmpigen:
	cd ncgen && $(MAKE)

ncmpidump:
	cd ncdump && $(MAKE)

clean:
	set -e; for i in $(SUBDIRS); do ( cd $$i && $(MAKE) clean ; ) ; done

install: ncgen_install ncdump_install

ncgen_install:
	cd ncgen && $(MAKE) install

ncdump_install:
	cd ncdump && $(MAKE) install
