if DEBUG
AM_CXXFLAGS=-O0 -g
else
AM_CXXFLAGS=-O3
endif

pkginclude_HEADERS=quickbam/bam.h\
				   quickbam/mbgzf.h\
				   quickbam/mfile.h\
				   quickbam/nfo_iterator.h\
				   quickbam/index.h\
				   quickbam/mpileup.h

lib_LTLIBRARIES=libquickbam.la

libquickbam_la_SOURCES = mbgzf.cc bam.cc index.cc
libquickbam_la_LDFLAGS = -version-info: 0.0.1
