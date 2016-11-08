AC_INIT([fastqcheck], m4_esyscmd([tr -d '\n' < VERSION]), [rd@sanger.ac.uk], [$(PACKAGE_NAME)$(AC_PACKAGE_VERSION)], [https://github.com/VertebrateResequencing/fastqcheck])

AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])   

AC_PROG_LIBTOOL
AC_PROG_CC                                   
AC_PROG_CXX

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

AC_CHECK_HEADERS([math.h])
