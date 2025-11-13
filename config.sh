#!/bin/bash

rm -rf CMakeFiles *Cache.txt *.cmake

COMPILER_PATH='/usr/bin/'

cmake \
    -DCMAKE_CXX_COMPILER:PATH=$COMPILER_PATH/g++ \
    -DCMAKE_C_COMPILER:PATH=$COMPILER_PATH/gcc \
    -DCMAKE_CXX_FLAGS:STRING='-Wall -fPIC -pipe' \
    -DCMAKE_C_FLAGS:STRING='-fPIC -pipe -Wall' \
    -DTPL_DIR:STRING="$HOME/tpl"

make -j16
