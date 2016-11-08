# ship these headers in the "make dist" target
dist_data_DATA = readseq.h

# fastqcheck is our top level program
bin_PROGRAMS = fastqcheck
fastqcheck_SOURCES = fastqcheck.c
fastqcheck_LDADD=-lm libfastqcheck.la

# libfastqcheck.so is our library
lib_LTLIBRARIES = libfastqcheck.la
libfastqcheck_la_SOURCES = readseq.c
libfastqcheck_la_LDFLAGS= -version-info 1:0
