.. The API section of the libquickbam documentation

Overview
========

This section of the document describe the types and interfaces of libquickbam
in details. It is meant to be a programming reference, instead of a proper
introduction on how to use libquickbam. For such an introduction, please see
Introduction.

The libquickbam API consists of several C structs and functions that either
take these structs as arguments, or return them as results. Some are meant to
be used as an opaque data type, whose internals are only meaningful to offered
library functions; while others are deliberately designed to mirror the byte
layout of data structures described in the SAM file format specification, and
meant to be accessed directly by client code. Here we offer a general overview
of these structs and their purposes.

Low level structs and functions
===============================

nfo_iterator
^^^^^^^^^^^^

.. doxygenstruct:: nfo_iterator
   :members:

file_slicer_t
^^^^^^^^^^^^^

.. doxygenstruct:: file_slicer_t
   :members:


bgzf_block_t
^^^^^^^^^^^^

Per SAM file specification, a BAM file consists of many BGZF compression
blocks. This struct is meant to mirror the byte layout of these blocks.
However, since BGZF blocks contain variable length fields, the struct can only
provide direct access to the fixed length fields, as well as the beginning of
the variable length field. It further made the assumption that each BGZF block
only contains a single "Extra" field with (66, 67) identifiers.

This struct is used by other components of libquickbam, and is rarely used in
client programs.

.. doxygenstruct:: bgzf_block_t
.. doxygenfunction:: is_bgzf_eof_block
.. doxygentypedef:: bgzf_iterator_t
.. doxygenfunction:: bgzf_isize
.. doxygenfunction:: bgzf_inflate
.. doxygenfunction:: bgzf_inflate_range
.. doxygenfunction:: bgzf_inflate_range_p

bgzf_slicer_iterator_t
^^^^^^^^^^^^^^^^^^^^^^

.. doxygenstruct:: bgzf_slicer_iterator_t
   :members:


BAM file related structs and functions
======================================

bam_header_t
^^^^^^^^^^^^

This struct mirrors the byte layout of the decompressed BAM file header, and
provides direct access to the fixed length fields. 

.. doxygenstruct:: bam_header_t
   :members:

.. doxygenfunction:: bam_buffer_contains_header

bam_rec_t
^^^^^^^^^

This struct mirrors the byte layout of the decompressed BAM file records, aka
sequence reads. It provides direct access to the fixed length fields. Functions
with ``bam_`` offers access to the fields encoded into the variable length
field.

.. doxygenstruct:: bam_rec_t
   :members:

.. doxygenfunction:: bam_load_block
.. doxygenfunction:: bam_count_records
.. doxygenfunction:: bam_query_length
.. doxygenfunction:: bam_read_name
.. doxygenfunction:: bam_cigar_ptr
.. doxygenfunction:: bam_seq_ptr
.. doxygenfunction:: bam_bqual_ptr
.. doxygenfunction:: bam_unpack_base
.. doxygenfunction:: bam_load_region


BAM index related structs and functions
=======================================

index_t
^^^^^^^^^^^

This struct mirrors the byte layout of the BAM index file format. This format
contains several layers of nested, variable length fields with their own byte
layouts, and are accessed by helper structs including ``index_ref_t``,
``index_bin_t``, and ``index_chunk_t``. Refer to the detailed API documentation
section regarding these.

.. doxygenstruct:: index_t
   :members:
.. doxygenstruct:: index_ref_t
   :members:
.. doxygenstruct:: index_bin_t
   :members:
.. doxygenstruct:: index_chunk_t
   :members:

.. doxygenfunction:: index_coffset
.. doxygenfunction:: index_uoffset
.. doxygenfunction:: index_read
.. doxygenfunction:: index_free
.. doxygenfunction:: index_to_regions


Multiple input pileup engine
============================

mpileup_t
^^^^^^^^^

This strcut is used by the mpileup engine to describe the pileup information at
a given genomic position. It is passed as the argument to the visitor function
provided by the client program when the pileup engine is called. As part of
this struct, another struct ``pileup_info_t`` is used to contain per-read
information at the piled up position.

.. doxygenstruct:: mpileup_t
   :members:

.. doxygenstruct:: pileup_info_t
   :members:

.. doxygentypedef:: indices_t
.. doxygenfunction:: mpileup
