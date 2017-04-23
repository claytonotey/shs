AC_INIT([shs],[0.0.0])
dnl Use the m4/ directory to contain libtool macros that will be needed to 
dnl go with the ltmain.sh script (both will be created/updated via the command
dnl libtoolize --copy --force
AC_CONFIG_MACRO_DIR([m4])

#------------------------------------------------------------------------------------
# Library's "official" release version:

LIBSHS_RELEASE_MAJOR=0
LIBSHS_RELEASE_MINOR=0
LIBSHS_RELEASE_BUILD=0
AC_CONFIG_SRCDIR([src/shs.cpp])
AC_CANONICAL_TARGET([])
AM_CONFIG_HEADER(src/config.h)

AC_DISABLE_SHARED dnl allows for optimizations

#------------------------------------------------------------------------------------
# The following is the libtool / shared library version. This doesn't have to
# do anything with the release version. It MUST conform to the following rules:
#
#  1. Start with version information of `0:0:0' for each libtool library.
#  2. Update the version information only immediately before a public release of
#     your software. More frequent updates are unnecessary, and only guarantee
#     that the current interface number gets larger faster.
#  3. If the library source code has changed at all since the last update, then
#     increment revision (`c:r:a' becomes `c:r+1:a').
#  4. If any interfaces have been added, removed, or changed since the last update,
#     increment current, and set revision to 0.
#  5. If any interfaces have been added since the last public release, then increment
#     age.
#  6. If any interfaces have been removed since the last public release, then set age
#     to 0.

SHARED_VERSION_INFO="0:0:0"


AM_PROG_LIBTOOL

cxxflags_save="$CXXFLAGS"
AC_LANG([C++])
AC_PROG_CXX
AC_PROG_F77
AX_PROG_CXX_MPI
CXXFLAGS="$cxxlags_save"


AC_SUBST(SHLIB_VERSION_ARG)
AC_SUBST(SHARED_VERSION_INFO)

AC_ARG_ENABLE(static,[  --enable-static enable static build: (default: disable)],[enable_static=$enableval],[enable_static=no])

AC_ARG_ENABLE(debug,[  --enable-debug enable debug build: (default: disable)],[enable_debug=$enableval],[enable_debug=no])

AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
         [compile with MPI (parallelization) support. If none is found,
             MPI is not used. Default: auto])
        ],,[with_mpi=auto])

AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[have_mpi=yes],[
        have_mpi=no
        if test x"$with_mpi" = xyes; then
          AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
        else
          AC_MSG_WARN([No MPI compiler found, won't use MPI.])
        fi
    ])

AC_ARG_WITH(omp, [AS_HELP_STRING([--with-omp],
            [compile with OpenMP support. If none is found,
                OpenMP is not used. Default: auto])
        ],,[with_omp=auto])

AX_OPENMP([have_omp=yes],[
        have_omp=no
        if test x"$with_omp" = xyes; then
          AC_MSG_FAILURE([OMP requested, but couldn't use OMP.])
        else
          AC_MSG_WARN([No OMP found, won't use OMP.])
        fi
    ])

SHS_CFLAGS="-I./util -I../"

if test x$enable_static = xno; then
  case "$target_os" in
     darwin*)
	SHS_CFLAGS="$SHS_CFLAGS -mdynamic-no-pic"
	;;
	*)
	;;
   esac
fi

if test x$enable_debug = xyes; then
   SHS_CFLAGS="$SHS_CFLAGS -g -O3"
else
   SHS_CFLAGS="$SHS_CFLAGS -O3 -fomit-frame-pointer"
fi


if test x$have_omp = xyes; then
   SHS_CFLAGS="$SHS_CFLAGS -fopenmp"
   AC_DEFINE(HAVE_OMP,1,[Define to use OpenMP])
fi

if test x$have_mpi = xyes; then
   AC_DEFINE(HAVE_MPI,1,[Define to use MPI])
fi

AM_CONDITIONAL(HAVE_MPI, test x$have_mpi = xyes)

AC_SUBST(SHS_CFLAGS)

AM_INIT_AUTOMAKE([])

AM_MAINTAINER_MODE()

AC_CONFIG_FILES([Makefile src/Makefile src/amos/Makefile src/util/Makefile libshs.spec shs.pc])
AC_OUTPUT
