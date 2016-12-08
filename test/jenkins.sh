#!/bin/sh
#
# MCS hosts a jenkins server (a continuous-integration build bot) at
# https://jenkins-ci.mcs.anl.gov/job/Parallel-NetCDF/

# this is just a beginning.
# - Would like to integrate coverity into the reports
# - our tests should report something Jenkins understands so that we get back
#   a finer-grained "N of M tests failed" instead of a single "pass/fail" for
#   all tests

# presumes the presence of the autotools and an mpi compiler.  The jenkins
# build slaves use the MCS workstation environment, which today is
# autoconf-2.68 and mpich-1.4.1
autoreconf -fi

# We are using the following Jenkins environment variables:
# WORKSPACE
#     The absolute path of the directory assigned to the build as a workspace.
# BUILD_TAG
#    String of "jenkins-${JOB_NAME}-${BUILD_NUMBER}". Convenient to put into a
#    resource file, a jar file, etc for easier identification.
./configure --prefix=${WORKSPACE:-`pwd`}/install-${BUILD_TAG:-`date +"Y%m%d-%H%M%S"`}
make clean
make
make testing
make ptest
