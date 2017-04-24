SUBDIRS = amos util

AM_CXXFLAGS = @SHS_CFLAGS@

libshsincludedir = $(includedir)

lib_LTLIBRARIES = libshs.la
libshs_la_SOURCES = SparseBLAS.cpp Bessel.cpp zpmn.cpp 

shs_LDADD = libshs.la ./util/libutil.la amos/libamos.a -lblas -llapack
shs_LDFLAGS = $(os_ldflags)

shsmpi_LDADD = libshs.la ./util/libutil.la amos/libamos.la -lblas -llapack
shsmpi_LDFLAGS = $(os_ldflags)

shs_SOURCES = SHS.cpp SphereHalfspace.cpp
shsmpi_SOURCES = SHSMPI.cpp SphereHalfspaceMPI.cpp

if HAVE_MPI
bin_PROGRAMS = shs shsmpi
else
bin_PROGRAMS = shs
endif