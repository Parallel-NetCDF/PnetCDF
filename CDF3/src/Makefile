all:
	cd lib && $(MAKE)
	if [ yes = yes ] ; then cd libf && $(MAKE) ; fi
	cd utils && $(MAKE)

distclean: clean
clean:
	cd lib && $(MAKE) clean
	cd libf && $(MAKE) clean
	cd utils && $(MAKE) clean

install:
	cd lib && $(MAKE) install
	if [ yes = yes ] ; then cd libf && $(MAKE) install ; fi
	cd utils && $(MAKE) install
