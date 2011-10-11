# Makefile.in generated by automake 1.11.1 from Makefile.am.
# src/python/rf/Makefile.  Generated from Makefile.in by configure.

# Copyright (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
# 2003, 2004, 2005, 2006, 2007, 2008, 2009  Free Software Foundation,
# Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.



#
# Copyright (C) 2009  The University of Texas at Austin.
# Copyright 2004 Free Software Foundation, Inc.
# 
# This file is part of Hydra: A wireless multihop testbed.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# 

# -*- Makefile -*-
#
# Copyright (C) 2009  The University of Texas at Austin.
# Copyright 2004,2006 Free Software Foundation, Inc.
# 
# This file is part of Hydra: A wireless multihop testbed.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# 

pkgdatadir = $(datadir)/gr-hydra
pkgincludedir = $(includedir)/gr-hydra
pkglibdir = $(libdir)/gr-hydra
pkglibexecdir = $(libexecdir)/gr-hydra
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = i686-pc-linux-gnu
host_triplet = i686-pc-linux-gnu
target_triplet = i686-pc-linux-gnu
DIST_COMMON = $(noinst_PYTHON) $(rfpython_PYTHON) \
	$(srcdir)/Makefile.am $(srcdir)/Makefile.in \
	$(srcdir)/run_tests.in $(top_srcdir)/Makefile.common
subdir = src/python/rf
ACLOCAL_M4 = $(top_srcdir)/aclocal.m4
am__aclocal_m4_deps = $(top_srcdir)/config/acx_pthread.m4 \
	$(top_srcdir)/config/gr_boost.m4 \
	$(top_srcdir)/config/gr_gprof.m4 \
	$(top_srcdir)/config/gr_no_undefined.m4 \
	$(top_srcdir)/config/gr_omnithread.m4 \
	$(top_srcdir)/config/gr_pwin32.m4 \
	$(top_srcdir)/config/gr_python.m4 \
	$(top_srcdir)/config/gr_scripting.m4 \
	$(top_srcdir)/config/gr_swig.m4 \
	$(top_srcdir)/config/gr_x86_64.m4 \
	$(top_srcdir)/config/grc_hydra_core.m4 \
	$(top_srcdir)/config/grc_hydra_gui.m4 \
	$(top_srcdir)/config/grc_hydra_mac.m4 \
	$(top_srcdir)/config/grc_hydra_mpif.m4 \
	$(top_srcdir)/config/grc_hydra_patch.m4 \
	$(top_srcdir)/config/grc_hydra_phy.m4 \
	$(top_srcdir)/config/grc_hydra_rf.m4 \
	$(top_srcdir)/config/grc_hydra_util.m4 \
	$(top_srcdir)/config/lf_cxx.m4 \
	$(top_srcdir)/config/lf_warnings.m4 \
	$(top_srcdir)/config/pkg.m4 $(top_srcdir)/configure.ac
am__configure_deps = $(am__aclocal_m4_deps) $(CONFIGURE_DEPENDENCIES) \
	$(ACLOCAL_M4)
mkinstalldirs = $(install_sh) -d
CONFIG_HEADER = $(top_builddir)/config.h
CONFIG_CLEAN_FILES = run_tests
CONFIG_CLEAN_VPATH_FILES =
SOURCES =
DIST_SOURCES =
am__vpath_adj_setup = srcdirstrip=`echo "$(srcdir)" | sed 's|.|.|g'`;
am__vpath_adj = case $$p in \
    $(srcdir)/*) f=`echo "$$p" | sed "s|^$$srcdirstrip/||"`;; \
    *) f=$$p;; \
  esac;
am__strip_dir = f=`echo $$p | sed -e 's|^.*/||'`;
am__install_max = 40
am__nobase_strip_setup = \
  srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*|]/\\\\&/g'`
am__nobase_strip = \
  for p in $$list; do echo "$$p"; done | sed -e "s|$$srcdirstrip/||"
am__nobase_list = $(am__nobase_strip_setup); \
  for p in $$list; do echo "$$p $$p"; done | \
  sed "s| $$srcdirstrip/| |;"' / .*\//!s/ .*/ ./; s,\( .*\)/[^/]*$$,\1,' | \
  $(AWK) 'BEGIN { files["."] = "" } { files[$$2] = files[$$2] " " $$1; \
    if (++n[$$2] == $(am__install_max)) \
      { print $$2, files[$$2]; n[$$2] = 0; files[$$2] = "" } } \
    END { for (dir in files) print dir, files[dir] }'
am__base_list = \
  sed '$$!N;$$!N;$$!N;$$!N;$$!N;$$!N;$$!N;s/\n/ /g' | \
  sed '$$!N;$$!N;$$!N;$$!N;s/\n/ /g'
am__installdirs = "$(DESTDIR)$(rfpythondir)"
py_compile = $(top_srcdir)/py-compile
am__tty_colors = \
red=; grn=; lgn=; blu=; std=
DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
ACLOCAL = ${SHELL} /home/testbed/hydra-0.4/gr-hydra/missing --run aclocal-1.11
AMTAR = ${SHELL} /home/testbed/hydra-0.4/gr-hydra/missing --run tar
AR = ar
AS = as
AUTOCONF = ${SHELL} /home/testbed/hydra-0.4/gr-hydra/missing --run autoconf
AUTOHEADER = ${SHELL} /home/testbed/hydra-0.4/gr-hydra/missing --run autoheader
AUTOMAKE = ${SHELL} /home/testbed/hydra-0.4/gr-hydra/missing --run automake-1.11
AWK = gawk
BOOST_CFLAGS = 
CC = gcc
CCAS = gcc
CCASDEPMODE = depmode=gcc3
CCASFLAGS = -g -O2
CCDEPMODE = depmode=gcc3
CFLAGS = -g -O2 -pthread
CPP = gcc -E
CPPFLAGS = 
CXX = g++
CXXCPP = g++ -E
CXXDEPMODE = depmode=gcc3
CXXFLAGS = -g -O2 -Wall -Woverloaded-virtual -pthread
CXX_FOR_BUILD = g++
CYGPATH_W = echo
DEFINES = 
DEFS = -DHAVE_CONFIG_H
DEPDIR = .deps
DLLTOOL = dlltool
DSYMUTIL = 
DUMPBIN = 
ECHO_C = 
ECHO_N = -n
ECHO_T = 
EGREP = /bin/grep -E
EXEEXT = 
FGREP = /bin/grep -F
GNURADIOSRCDIR = /home/testbed/hydra-0.4/gnuradio-3.2.2
GNURADIO_CORE_CFLAGS = -pthread -DOMNITHREAD_POSIX=1 -I/home/testbed/hydra-0.4/gr/include/gnuradio -I/home/testbed/hydra-0.4/gr/include  
GNURADIO_CORE_INCLUDEDIR = /home/testbed/hydra-0.4/gr/include/gnuradio
GNURADIO_CORE_LIBS = -L/home/testbed/hydra-0.4/gr/lib -lgnuradio-core -lboost_thread-mt -lrt -lboost_date_time-mt -lgruel -lfftw3f -lgsl -lgslcblas -lm -lgromnithread  
GREP = /bin/grep
GUILE = /usr/bin/guile
INSTALL = /usr/bin/install -c
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = $(install_sh) -c -s
ITPP_CFLAGS = -Wall -g -DNDEBUG -pipe -O2  
ITPP_INCLUDEDIR = /usr/include
ITPP_LIBS = -litpp  
LD = /usr/bin/ld
LDFLAGS = 
LIBOBJS = 
LIBS =  -L/home/testbed/hydra-0.4/gr/lib -lgnuradio-core -lboost_thread-mt -lrt -lboost_date_time-mt -lgruel -lfftw3f -lgsl -lgslcblas -lm -lgromnithread  
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LIPO = 
LN_S = ln -s
LTLIBOBJS = 
MAKEINFO = ${SHELL} /home/testbed/hydra-0.4/gr-hydra/missing --run makeinfo
MKDIR_P = /bin/mkdir -p
NM = /usr/bin/nm -B
NMEDIT = 
NO_UNDEFINED = 
OBJDUMP = objdump
OBJEXT = o
OTOOL = 
OTOOL64 = 
PACKAGE = gr-hydra
PACKAGE_BUGREPORT = 
PACKAGE_NAME = 
PACKAGE_STRING = 
PACKAGE_TARNAME = 
PACKAGE_URL = 
PACKAGE_VERSION = 
PATH_SEPARATOR = :
PKG_CONFIG = /usr/bin/pkg-config
PTHREAD_CC = gcc
PTHREAD_CFLAGS = -pthread
PTHREAD_LIBS =  -lrt 
PYTHON = /usr/bin/python
PYTHON_CPPFLAGS = -I/usr/include/python2.6
PYTHON_EXEC_PREFIX = ${exec_prefix}
PYTHON_LDFLAGS = 
PYTHON_PLATFORM = linux2
PYTHON_PREFIX = ${prefix}
PYTHON_VERSION = 2.6
RANLIB = ranlib
RM_PROG = /bin/rm
SED = /bin/sed
SET_MAKE = 
SHELL = /bin/sh
STD_DEFINES_AND_INCLUDES = -pthread -DOMNITHREAD_POSIX=1 -I/home/testbed/hydra-0.4/gr/include/gnuradio -I/home/testbed/hydra-0.4/gr/include   
STRIP = strip
SWIG = /usr/bin/swig -c++
SWIG_PYTHON_CPPFLAGS = -I/usr/include/python2.6
SWIG_PYTHON_LIB = -lswigpy
SWIG_PYTHON_OPT = -python
VERSION = 0.1svn
abs_builddir = /home/testbed/hydra-0.4/gr-hydra/src/python/rf
abs_srcdir = /home/testbed/hydra-0.4/gr-hydra/src/python/rf
abs_top_builddir = /home/testbed/hydra-0.4/gr-hydra
abs_top_srcdir = /home/testbed/hydra-0.4/gr-hydra
ac_ct_CC = gcc
ac_ct_CXX = g++
ac_ct_DUMPBIN = 
am__include = include
am__leading_dot = .
am__quote = 
am__tar = ${AMTAR} chof - "$$tardir"
am__untar = ${AMTAR} xf -
bindir = ${exec_prefix}/bin
build = i686-pc-linux-gnu
build_alias = 
build_cpu = i686
build_os = linux-gnu
build_vendor = pc
builddir = .
datadir = ${datarootdir}
datarootdir = ${prefix}/share
docdir = ${datarootdir}/doc/${PACKAGE}
dvidir = ${docdir}
exec_prefix = ${prefix}
host = i686-pc-linux-gnu
host_alias = 
host_cpu = i686
host_os = linux-gnu
host_vendor = pc
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = ${SHELL} /home/testbed/hydra-0.4/gr-hydra/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
lt_ECHO = echo
mandir = ${datarootdir}/man
mkdir_p = /bin/mkdir -p
oldincludedir = /usr/include
pdfdir = ${docdir}
pkgpyexecdir = ${pyexecdir}/gr-hydra
pkgpythondir = ${pythondir}/gr-hydra
prefix = /home/testbed/hydra-0.4/gr
program_transform_name = s,x,x,
psdir = ${docdir}
pyexecdir = ${exec_prefix}/lib/python2.6/site-packages
pythondir = ${prefix}/lib/python2.6/site-packages
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
srcdir = .
sysconfdir = ${prefix}/etc
target = i686-pc-linux-gnu
target_alias = 
target_cpu = i686
target_os = linux-gnu
target_vendor = pc
top_build_prefix = ../../../
top_builddir = ../../..
top_srcdir = ../../..

# usual gnuradio definitions
grincludedir = $(includedir)/gnuradio
grswigincludedir = $(grincludedir)/swig
grpythondir = $(pythondir)/gnuradio
grpyexecdir = $(pyexecdir)/gnuradio

# swig flags
SWIGPYTHONFLAGS = -fvirtual -python -modern
SWIGGRFLAGS = -I$(GNURADIO_CORE_INCLUDEDIR)/swig -I$(GNURADIO_CORE_INCLUDEDIR)

# Don't assume that make predefines $(RM), because BSD make does
# not. We define it now in configure.ac using AM_PATH_PROG, but now
# here have to add a -f to be like GNU make.
RM = $(RM_PROG) -f

# Additional definitions for hydra

# includes
# This ends up at:
# 	${prefix}/include/hydra[/swig]
hydraincludedir = $(includedir)/hydra
hydraswigincludedir = $(hydraincludedir)/swig

# python install path
# This ends up at:
#   ${prefix}/lib/python${python_version}/site-packages/hydra
hydrapythondir = $(pythondir)/hydra
hydrapyexecdir = $(pyexecdir)/hydra

# itpp flags
ITPPCFLAGS = `itpp-config --cflags`
ITPPLDFLAGS = `itpp-config --libs`

# included to build using guile
#RUN_GUILE = GUILE_LOAD_PATH="/home/testbed/hydra-0.4/gr-hydra/pmt/src/scheme:/home/testbed/hydra-0.4/gr-hydra/mblock/src/scheme" /usr/bin/guile -e main -s
RUN_GUILE = GUILE_LOAD_PATH="$(GNURADIOSRCDIR)/pmt/src/scheme:$(GNURADIOSRCDIR)/mblock/src/scheme" /usr/bin/guile -e main -s
#COMPILE_MBH = $(RUN_GUILE) $(top_srcdir)/mblock/src/scheme/gnuradio/compile-mbh.scm
COMPILE_MBH = $(RUN_GUILE) $(GNURADIOSRCDIR)/mblock/src/scheme/gnuradio/compile-mbh.scm
DISTCLEANFILES = Makefile.in
EXTRA_DIST = run_tests.in
rfpythondir = $(hydrapythondir)/PyHydra/rf
rfpython_PYTHON = \
	__init__.py	\
	radiorx.py	\
	radiotx.py

TESTS = \
	run_tests

noinst_PYTHON = 
CLEANFILES = *.pyc *.pyo
all: all-am

.SUFFIXES:
$(srcdir)/Makefile.in:  $(srcdir)/Makefile.am $(top_srcdir)/Makefile.common $(am__configure_deps)
	@for dep in $?; do \
	  case '$(am__configure_deps)' in \
	    *$$dep*) \
	      ( cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh ) \
	        && { if test -f $@; then exit 0; else break; fi; }; \
	      exit 1;; \
	  esac; \
	done; \
	echo ' cd $(top_srcdir) && $(AUTOMAKE) --gnu src/python/rf/Makefile'; \
	$(am__cd) $(top_srcdir) && \
	  $(AUTOMAKE) --gnu src/python/rf/Makefile
.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@case '$?' in \
	  *config.status*) \
	    cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh;; \
	  *) \
	    echo ' cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe)'; \
	    cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe);; \
	esac;

$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh

$(top_srcdir)/configure:  $(am__configure_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(ACLOCAL_M4):  $(am__aclocal_m4_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(am__aclocal_m4_deps):
run_tests: $(top_builddir)/config.status $(srcdir)/run_tests.in
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs
install-rfpythonPYTHON: $(rfpython_PYTHON)
	@$(NORMAL_INSTALL)
	test -z "$(rfpythondir)" || $(MKDIR_P) "$(DESTDIR)$(rfpythondir)"
	@list='$(rfpython_PYTHON)'; dlist=; list2=; test -n "$(rfpythondir)" || list=; \
	for p in $$list; do \
	  if test -f "$$p"; then b=; else b="$(srcdir)/"; fi; \
	  if test -f $$b$$p; then \
	    $(am__strip_dir) \
	    dlist="$$dlist $$f"; \
	    list2="$$list2 $$b$$p"; \
	  else :; fi; \
	done; \
	for file in $$list2; do echo $$file; done | $(am__base_list) | \
	while read files; do \
	  echo " $(INSTALL_DATA) $$files '$(DESTDIR)$(rfpythondir)'"; \
	  $(INSTALL_DATA) $$files "$(DESTDIR)$(rfpythondir)" || exit $$?; \
	done || exit $$?; \
	if test -n "$$dlist"; then \
	  if test -z "$(DESTDIR)"; then \
	    PYTHON=$(PYTHON) $(py_compile) --basedir "$(rfpythondir)" $$dlist; \
	  else \
	    PYTHON=$(PYTHON) $(py_compile) --destdir "$(DESTDIR)" --basedir "$(rfpythondir)" $$dlist; \
	  fi; \
	else :; fi

uninstall-rfpythonPYTHON:
	@$(NORMAL_UNINSTALL)
	@list='$(rfpython_PYTHON)'; test -n "$(rfpythondir)" || list=; \
	files=`for p in $$list; do echo $$p; done | sed -e 's|^.*/||'`; \
	test -n "$$files" || exit 0; \
	filesc=`echo "$$files" | sed 's|$$|c|'`; \
	fileso=`echo "$$files" | sed 's|$$|o|'`; \
	echo " ( cd '$(DESTDIR)$(rfpythondir)' && rm -f" $$files ")"; \
	cd "$(DESTDIR)$(rfpythondir)" && rm -f $$files || exit $$?; \
	echo " ( cd '$(DESTDIR)$(rfpythondir)' && rm -f" $$filesc ")"; \
	cd "$(DESTDIR)$(rfpythondir)" && rm -f $$filesc || exit $$?; \
	echo " ( cd '$(DESTDIR)$(rfpythondir)' && rm -f" $$fileso ")"; \
	cd "$(DESTDIR)$(rfpythondir)" && rm -f $$fileso
tags: TAGS
TAGS:

ctags: CTAGS
CTAGS:


check-TESTS: $(TESTS)
	@failed=0; all=0; xfail=0; xpass=0; skip=0; \
	srcdir=$(srcdir); export srcdir; \
	list=' $(TESTS) '; \
	$(am__tty_colors); \
	if test -n "$$list"; then \
	  for tst in $$list; do \
	    if test -f ./$$tst; then dir=./; \
	    elif test -f $$tst; then dir=; \
	    else dir="$(srcdir)/"; fi; \
	    if $(TESTS_ENVIRONMENT) $${dir}$$tst; then \
	      all=`expr $$all + 1`; \
	      case " $(XFAIL_TESTS) " in \
	      *[\ \	]$$tst[\ \	]*) \
		xpass=`expr $$xpass + 1`; \
		failed=`expr $$failed + 1`; \
		col=$$red; res=XPASS; \
	      ;; \
	      *) \
		col=$$grn; res=PASS; \
	      ;; \
	      esac; \
	    elif test $$? -ne 77; then \
	      all=`expr $$all + 1`; \
	      case " $(XFAIL_TESTS) " in \
	      *[\ \	]$$tst[\ \	]*) \
		xfail=`expr $$xfail + 1`; \
		col=$$lgn; res=XFAIL; \
	      ;; \
	      *) \
		failed=`expr $$failed + 1`; \
		col=$$red; res=FAIL; \
	      ;; \
	      esac; \
	    else \
	      skip=`expr $$skip + 1`; \
	      col=$$blu; res=SKIP; \
	    fi; \
	    echo "$${col}$$res$${std}: $$tst"; \
	  done; \
	  if test "$$all" -eq 1; then \
	    tests="test"; \
	    All=""; \
	  else \
	    tests="tests"; \
	    All="All "; \
	  fi; \
	  if test "$$failed" -eq 0; then \
	    if test "$$xfail" -eq 0; then \
	      banner="$$All$$all $$tests passed"; \
	    else \
	      if test "$$xfail" -eq 1; then failures=failure; else failures=failures; fi; \
	      banner="$$All$$all $$tests behaved as expected ($$xfail expected $$failures)"; \
	    fi; \
	  else \
	    if test "$$xpass" -eq 0; then \
	      banner="$$failed of $$all $$tests failed"; \
	    else \
	      if test "$$xpass" -eq 1; then passes=pass; else passes=passes; fi; \
	      banner="$$failed of $$all $$tests did not behave as expected ($$xpass unexpected $$passes)"; \
	    fi; \
	  fi; \
	  dashes="$$banner"; \
	  skipped=""; \
	  if test "$$skip" -ne 0; then \
	    if test "$$skip" -eq 1; then \
	      skipped="($$skip test was not run)"; \
	    else \
	      skipped="($$skip tests were not run)"; \
	    fi; \
	    test `echo "$$skipped" | wc -c` -le `echo "$$banner" | wc -c` || \
	      dashes="$$skipped"; \
	  fi; \
	  report=""; \
	  if test "$$failed" -ne 0 && test -n "$(PACKAGE_BUGREPORT)"; then \
	    report="Please report to $(PACKAGE_BUGREPORT)"; \
	    test `echo "$$report" | wc -c` -le `echo "$$banner" | wc -c` || \
	      dashes="$$report"; \
	  fi; \
	  dashes=`echo "$$dashes" | sed s/./=/g`; \
	  if test "$$failed" -eq 0; then \
	    echo "$$grn$$dashes"; \
	  else \
	    echo "$$red$$dashes"; \
	  fi; \
	  echo "$$banner"; \
	  test -z "$$skipped" || echo "$$skipped"; \
	  test -z "$$report" || echo "$$report"; \
	  echo "$$dashes$$std"; \
	  test "$$failed" -eq 0; \
	else :; fi

distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  sed -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   sed '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | sed -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done
check-am: all-am
	$(MAKE) $(AM_MAKEFLAGS) check-TESTS
check: check-am
all-am: Makefile
installdirs:
	for dir in "$(DESTDIR)$(rfpythondir)"; do \
	  test -z "$$dir" || $(MKDIR_P) "$$dir"; \
	done
install: install-am
install-exec: install-exec-am
install-data: install-data-am
uninstall: uninstall-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	  install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	  `test -z '$(STRIP)' || \
	    echo "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'"` install
mostlyclean-generic:

clean-generic:
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)

distclean-generic:
	-test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)
	-test . = "$(srcdir)" || test -z "$(CONFIG_CLEAN_VPATH_FILES)" || rm -f $(CONFIG_CLEAN_VPATH_FILES)
	-test -z "$(DISTCLEANFILES)" || rm -f $(DISTCLEANFILES)

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
clean: clean-am

clean-am: clean-generic clean-libtool mostlyclean-am

distclean: distclean-am
	-rm -f Makefile
distclean-am: clean-am distclean-generic

dvi: dvi-am

dvi-am:

html: html-am

html-am:

info: info-am

info-am:

install-data-am: install-rfpythonPYTHON

install-dvi: install-dvi-am

install-dvi-am:

install-exec-am:

install-html: install-html-am

install-html-am:

install-info: install-info-am

install-info-am:

install-man:

install-pdf: install-pdf-am

install-pdf-am:

install-ps: install-ps-am

install-ps-am:

installcheck-am:

maintainer-clean: maintainer-clean-am
	-rm -f Makefile
maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-am

mostlyclean-am: mostlyclean-generic mostlyclean-libtool

pdf: pdf-am

pdf-am:

ps: ps-am

ps-am:

uninstall-am: uninstall-rfpythonPYTHON

.MAKE: check-am install-am install-strip

.PHONY: all all-am check check-TESTS check-am clean clean-generic \
	clean-libtool distclean distclean-generic distclean-libtool \
	distdir dvi dvi-am html html-am info info-am install \
	install-am install-data install-data-am install-dvi \
	install-dvi-am install-exec install-exec-am install-html \
	install-html-am install-info install-info-am install-man \
	install-pdf install-pdf-am install-ps install-ps-am \
	install-rfpythonPYTHON install-strip installcheck \
	installcheck-am installdirs maintainer-clean \
	maintainer-clean-generic mostlyclean mostlyclean-generic \
	mostlyclean-libtool pdf pdf-am ps ps-am uninstall uninstall-am \
	uninstall-rfpythonPYTHON


# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT: