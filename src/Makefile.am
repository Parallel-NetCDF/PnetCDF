#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id$
#
# @configure_input@

# adding . to the list of SUBDIRS causes the directories built before .
# This is useful for "make -j" option, as libs should be built after all
# sub-libraries
SUBDIRS = include dispatchers drivers binding . libs utils packaging

dist-hook:
	-rm -rf `find $(distdir) -name .svn`

