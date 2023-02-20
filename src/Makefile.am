if DEBUG
AM_CXXFLAGS=-O0 -g
else
AM_CXXFLAGS=-O3
endif

pkginclude_HEADERS=mmbam/bam.h\
				   mmbam/mbgzf.h\
				   mmbam/mfile.h\
				   mmbam/nfo_iterator.h\
				   mmbam/index.h\
				   mmbam/mpileup.h

lib_LTLIBRARIES=libmmbam.la

libmmbam_la_SOURCES = mbgzf.cc bam.cc index.cc
libmmbam_la_LDFLAGS = -version-info: 0.0.1
