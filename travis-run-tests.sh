#!/bin/bash
set -ev

# Coverity Scan addons build_command_prepend already runs configure
# ./configure

# Coverity Scan addons build_command already runs make distcheck
# make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent"

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent --enable-shared"

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent --enable-shared --enable-burst-buffering"

