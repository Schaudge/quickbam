bin_PROGRAMS=readcount-mmbam flagstats-mmbam snp-pileup-mmbam

AM_CPPFLAGS = -I$(top_srcdir)/src

readcount_mmbam_SOURCES  = readcount-mmbam.cc
readcount_mmbam_CXXFLAGS = -fopenmp
readcount_mmbam_LDADD    = $(top_builddir)/src/libmmbam.la -ltbb -ldeflate

flagstats_mmbam_SOURCES = flagstats-mmbam.cc
flagstats_mmbam_LDADD   = $(top_builddir)/src/libmmbam.la -ltbb -ldeflate

snp_pileup_mmbam_SOURCES = snp-pileup-mmbam.cc
snp_pileup_mmbam_LDADD   = $(top_builddir)/src/libmmbam.la -ltbb -ldeflate
