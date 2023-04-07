if DEBUG
AM_CXXFLAGS=-O0 -g
else
AM_CXXFLAGS=-O3
endif

bin_PROGRAMS=readcount-quickbam flagstats-quickbam snp-pileup-quickbam

AM_CPPFLAGS = -I$(top_srcdir)/src

readcount_quickbam_SOURCES  = readcount-quickbam.cc
readcount_quickbam_CXXFLAGS = -fopenmp
readcount_quickbam_LDADD    = $(top_builddir)/src/libquickbam.la -ltbb -ldeflate

flagstats_quickbam_SOURCES = flagstats-quickbam.cc
flagstats_quickbam_CXXFLAGS = -fopenmp
flagstats_quickbam_LDADD   = $(top_builddir)/src/libquickbam.la -ltbb -ldeflate

snp_pileup_quickbam_SOURCES = snp-pileup-quickbam.cc
snp_pileup_quickbam_LDADD   = $(top_builddir)/src/libquickbam.la -ltbb -ltbbmalloc -ldeflate
