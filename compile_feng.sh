#! /bin/bash

module add qt/5.13.0
qmake Simplification.pro
make -j8

