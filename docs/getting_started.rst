API Structure at a high level
=============================

The entire libquickbam library API is designed with the following principal:

* use lightweight data structure as much as possible. Classes are mostly
  avoided

* For each data type ``type``, APIs that act on it are named ``<type>-<action>``

.. admonition:: Example

   function ``bam_buffer_contains_header(const
   std::vector<uint8_t>&)`` returns true if the buffer contains a complete BAM
   header starting from the beginning

A typical client code (programs that use libquickbam) follows these steps:

1. create an mfile_t object initialized with the BAM file the program intend to
   read

2. create an index_t object with the index to the BAM file to facilitate
   parallel reading

3. create parallel regions that can be processed independently, either from the
   index or other sources

.. admonition:: Example

   the ``readcount-quickbam`` and ``flagstats-quickbam`` example programs
   creates one region for every 16kb genome window, directly from the index

   the ``snp-pileup-quickbam`` example code creates one region for all variants
   whose locations span less than 1 megabyte of data in the BAM file, using the
   input VCF file as well as the BAM index

4. process the sequence reads in the BAM file in parallel across these regions

5. combine results from the parallel processing into final, global results

Compilation and Linking
=======================

libquickbam uses the GNU autotools build system. To compile and install libquickbam,
execute the following commands on the command line

.. code-block:: bash

  tar -xzf quickbam-<version>.tar.gz
  cd quickbam-<version>
  ./configure
  make
  sudo make install


libquickbam depends on libdeflate for decompression and Intel Thread Building
Block for parallelization. Programs that use libquickbam needs to specify
``-lquickbam -ldeflate -ltbb`` in their linker flags in order to properly link
these libraries.

Setting up a project
====================

We recommend using autotools or CMake to setup the build system of your new
project. If you prefer to write makefiles manually, make sure to 

* enable C++14 features (e.g. ``-std=c++14``) during compilation
* include the appropriate libraries ``-lquickbam -ldeflate -ltbb`` during linking.  

Example build systems using automake is available under the ``code_examples``
directory in the libquickbam repository.

Reading a BAM file without the index
====================================

Although the biggest advantage of using libquickbam, namely parallel reading,
cannot be achieved without an index file, it is still possible to iterate over
the BGZF blocks in a BAM file to e.g. parse the header. The following snippet
demonstrates opening a BAM file as a memory mapped file, and decompressing the
BGZF blocks until a complete BAM header is retrieved

Iterate over BGZF blocks using explicit loop
--------------------------------------------

.. code-block:: cpp
   
   #include <quickbam/mfile.h>
   #include <quickbam/mbgzf.h>
   #include <quickbam/bam.h>

   #include <string>

   int main(void) {
       std::string bam_filename="some/bam/file.bam";

       auto mfile = mfile_open(bam_filename);  // open BAM file as an mfile
       bgzf_mfile_proxy_t bgzf_proxy(mfile);   // create bgzf proxy over the mfile
       std::vector<uint8_t> buffer;            // decompression buffer

       for(auto& bgzf_block : bgzf_proxy) {    // iterate over all bgzf blocks
           auto block_buff = bgzf_inflate(bgzf_block); // decompress one block

           // append the decompressed bytes from this block to buffer
           buffer.insert(buffer.end(), block_buff.begin(), block_buff.end());

           // check if the current buffer contains a complete BAM header
           // if yes, break out of the loop
           if(bam_buffer_contains_header(buffer) break;
       }

       // no need to manually close the mfile since it is managed by a smart pointer
   }

Iterate over BGZF blocks using std::find_if
-------------------------------------------

In the next example, we achieve the same thing but with c++ std::find_if
instead of a explicit loop

.. code-block:: cpp
   
   // ... setup code ... //

   std::string bam_filename="some/bam/file.bam";

   auto mfile = mfile_open(bam_filename);  // open BAM file as an mfile
   bgzf_mfile_proxy_t bgzf_proxy(mfile);   // create bgzf proxy over the mfile
   std::vector<uint8_t> buffer;            // decompression buffer

   auto header_end = std::find_if(
       bgzf_proxy.begin(),
       bgzf_proxy.end(),
       [&buffer](auto& bgzf_block) {
           if(bam_buffer_contains_header(buffer)) return true;
           auto block_buffer = bgzf_inflate(bgzf_block);
           buffer.insert(buffer.end(), block_buffer.begin(), block_buffer.end());
           return false;
       }

    // at this point, header_end points at the BGZF block after the end of
    // header block, which will also be the beginning of the first bam record

    // ... rest of code ... //

Iterating and processing BAM records
====================================

Parsing a BAM record
--------------------

Once a byte vector of decompressed BAM records is acquired, we can use
bam_rec_t to parse the different data fields of a BAM record

.. code-block:: cpp

   // ... setup code ... //

   std::string bam_filename="some/bam/file.bam";

   auto mfile = mfile_open(bam_filename);  // open BAM file as an mfile
   bgzf_mfile_proxy_t bgzf_proxy(mfile);   // create bgzf proxy over the mfile
   std::vector<uint8_t> buffer;            // decompression buffer

   auto header_end = std::find_if(
       bgzf_proxy.begin(),
       bgzf_proxy.end(),
       [&buffer](auto& bgzf_block) {
           if(bam_buffer_contains_header(buffer)) return true;
           auto block_buffer = bgzf_inflate(bgzf_block);
           buffer.insert(buffer.end(), block_buffer.begin(), block_buffer.end());
           return false;
       }

    // at this point, header_end points at the BGZF block after the end of
    // header block, which will also be the beginning of the first bam record

    auto first_bam_block = bgzf_inflate(*header_end); // decompress first bam record block
    bam_rec_t* record = reinterpret_cast<bam_rec_t*>(first_bam_block.data());

    // print read name
    std::cout<<"READ NAME: "<<bam_read_name(record)<<std::endl;

    // find where the next read is
    record = BAM_NEXT(record);

    // you should check if record is still within buffer
    if(record < first_bam_block.data() + first_bam_block.size())
        std::cout<<"NEXT READ NAME: "<<bam_read_name(record)<<std::endl;

Iterating over a BAM records buffer
-----------------------------------

bam_iterator is designed to offer the ability to iterate over a buffer with
multiple BAM records. 

.. code-block:: cpp

   // ... setup code ...//
   // let's say buffer now contains multiple bam records

   // bam_iterator is a specialization of nfo_iterator_t
   // which can be initialized with a byte vector
   bam_iterator bam_it(buffer);
   bam_iterator bam_end(buffer, buffer.size());

   while(bam_it < bam_end) {
       // bam_iterator can be dereferenced and implicitly cast to bam_rec_t*
       std::cout<<"READ POS: "<<bam_it->ref_id<<":"<<bam_it->pos<<std::endl;
       std::cout<<"READ NAME: "<<bam_read_name(bam_it)<<std::endl;
       bam_it++;
   }

Random BAM accessing using the index
====================================

The BAM index file contains the necessary information to load sequence reads of
a given genomic region. Libquickbam takes advantage of the "linear index", which
contains the file offsets (compressed) and buffer offsets (decompressed) for
each 16kb genomic window. The following example demonstrates loading the BAM
records of a particular region

.. code-block:: cpp

   #include <quickbam/mfile.h>
   #include <quickbam/mbgzf.h>
   #include <quickbam/bam.h>
   #include <quickbam/index.h>

   #include <string>
   #include <fstream>

   int main(void) {
       std::string bam_filename="some/bam/file.bam";
       std::string bai_filename="some/bam/file.bam.bai";

       auto mfile = mfile_open(bam_filename);  // open BAM file as an mfile

       // open and parse index file
       auto bai_stream = std::ifstream(bai_filename);
       auto index = index_read(bai_stream);

       // reference contig name to ref_id map is in the BAM file header
       // for this example, we are hard coding ref_id to be 9
       // which most likely will correspond to chr10

       uint32_t region_start = 1500000;
       uint32_t region_end   = 1530000;

       auto buffer = bam_load_region(mfile, index, 9, region_start, region_end);

       // buffer now contains all reads on chromosome 10, between 1,500,000 bp
       // and 1,530,000 bp.

       bam_iterator bam_it(buffer);
       bam_iterator bam_end(buffer, buffer.size());

       while(bam_it < bam_end) {
           // process the read

           // advance the iterator
           bam_it++;
       }
   }
