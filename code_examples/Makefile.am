if DEBUG
AM_CXXFLAGS=-O0 -g
else
AM_CXXFLAGS=-O3
endif

bin_PROGRAMS=flagstats-quickbam snp-pileup-quickbam

if HAS_OPENMP
bin_PROGRAMS += readcount-quickbam 
endif

AM_CPPFLAGS = -I$(top_srcdir)/src

if HAS_OPENMP
readcount_quickbam_SOURCES  = readcount-quickbam.cc
readcount_quickbam_CXXFLAGS = @OPENMP_FLAGS@
readcount_quickbam_LDADD    = $(top_builddir)/src/libquickbam.la -ltbb -ldeflate
endif

flagstats_quickbam_SOURCES = flagstats-quickbam.cc
flagstats_quickbam_LDADD   = $(top_builddir)/src/libquickbam.la -ltbb -ldeflate

snp_pileup_quickbam_SOURCES = snp-pileup-quickbam.cc
snp_pileup_quickbam_LDADD   = $(top_builddir)/src/libquickbam.la -ltbb -ltbbmalloc -ldeflate
