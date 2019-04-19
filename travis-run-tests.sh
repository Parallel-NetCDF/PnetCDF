#!/bin/bash
set -ev

# Coverity Scan addons build_command_prepend already runs configure
# ./configure

# Coverity Scan addons build_command already runs make distcheck
# make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent"

# make -s V=1 LIBTOOLFLAGS=--silent -j4 tests

# make -s V=1 LIBTOOLFLAGS=--silent ptest

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent --enable-shared --enable-burst-buffering"

