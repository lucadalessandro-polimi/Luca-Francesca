#!/bin/bash

cpp=$HOME/fdaPDE-cpp
core=$cpp/fdaPDE/core
eigen=/usr/include/eigen3
include=include/
exe=$(basename -- $1 .cpp)

if [ -f $exe ]; then rm $exe; fi
g++ -o $exe $1 -I$cpp -I$core -I$eigen -I$include -O2 -std=c++20 -march=native -s
