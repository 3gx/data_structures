#!/bin/sh
gcc -c rtc.c
# ifort    -O4 -o polyhdb polyhdb.f -O4 -vec-report rtc.o -cpp
gfortran-mp-4.5 -O3 -o polyhdb polyhdb.f -O4 rtc.o -cpp
./polyhdb < rand.in


