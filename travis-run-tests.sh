#!/bin/bash
set -ev

./configure

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent"

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent --enable-shared"

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent --enable-shared --disable-static"

