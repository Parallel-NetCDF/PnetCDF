#!/bin/bash
set -ev

./configure

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS=--silent

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent --enable-shared LDFLAGS='-Wl,--allow-shlib-undefined'"

make distcheck -s V=1 LIBTOOLFLAGS=--silent DISTCHECK_CONFIGURE_FLAGS="--silent --enable-shared --disable-static LDFLAGS='-Wl,--allow-shlib-undefined'"

