dnl If anyone wants to improve on this go right ahead
AC_INIT([xAtlas], 0.1)
AC_PROG_CXX
AC_PREREQ([2.63])
AC_CONFIG_HEADERS([config.h])

dnl Where is htslib
AC_ARG_WITH([htslib_include],
            [AC_HELP_STRING([--with-htslib-include=DIR],
                            [HTSlib include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])

AC_ARG_WITH([htslib_lib],
            [AC_HELP_STRING([--with-htslib-lib=DIR],
                            [HTSlib library directory])],
            [LDFLAGS="-L$withval $CPPFLAGS"])

AC_ARG_WITH([htslib],
            [AC_HELP_STRING([--with-htslib=DIR],
                            [HTSlib installation directory])],
            [
            if test -z "$with_htslib_lib" -a -z "$with_htslib_include" ; then
                CPPFLAGS="-I$withval/include $CPPFLAGS"
                LDFLAGS="-L$withval/lib $LDFLAGS"
            else
                AC_MSG_FAILURE([Do not use --with-htslib and --with-htslib-include/--with-htslib-lib options simultaneously.])
            fi
            ])

dnl Where is pthread
AC_ARG_WITH([pthread_include],
            [AC_HELP_STRING([--with-pthread-include=DIR],
                            [pthread include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])

AC_ARG_WITH([pthread_lib],
            [AC_HELP_STRING([--with-pthread-lib=DIR],
                            [pthread library directory])],
            [LDFLAGS="-L$withval $CPPFLAGS"])

AC_ARG_WITH([pthread],
            [AC_HELP_STRING([--with-pthread=DIR],
                            [pthread installation directory])],
            [
            if test -z "$with_pthread_lib" -a -z "$with_pthread_include" ; then
                CPPFLAGS="-I$withval/include $CPPFLAGS"
                LDFLAGS="-L$withval/lib $LDFLAGS"
            else
                AC_MSG_FAILURE([Do not use --with-pthread and --with-pthread-include/--with-pthread-lib options simultaneously.])
            fi
            ])

dnl Check zlib? htslib should already require this
dnl AC_CHECK_HEADER([zlib.h],
dnl     [AC_SEARCH_LIBS([deflate], [z],
dnl         [],
dnl         [AC_MSG_ERROR([Failed to find zlib library])])],
dnl     [AC_MSG_ERROR([Failed to find zlib headers])])

dnl Check htslib
AC_CHECK_HEADER([htslib/hts.h],
    [AC_SEARCH_LIBS([hts_open], [hts],
        [],
        [AC_MSG_ERROR([Failed to find HTSLib library])])],
    [AC_MSG_ERROR([Failed to find HTSLib headers])])

dnl Check pthread/pthread_barrier extension needed for -p option
AC_ARG_ENABLE([multithreading],
    [AS_HELP_STRING([--enable-multithreading],
                    [enable multithreading (requires pthread_barrier)])])

AS_IF([test "x$enable_multithreading" = "xyes"], [
    AC_CHECK_HEADER([pthread.h],
        [AC_SEARCH_LIBS([pthread_create], [pthread],
            [AC_SEARCH_LIBS([pthread_barrier_init], [pthread],
                [AC_DEFINE([USE_PTHREAD])],
                [AC_MSG_WARN([Failed to find function pthread_barrier_init, disabling xAtlas multithreading support])])],
            [AC_MSG_WARN([Failed to find pthread library, disabling xAtlas multithreading support])])],
        [AC_MSG_WARN([Failed to find pthread headers, disabling xAtlas multithreading support])])
])

AC_OUTPUT([Makefile])
