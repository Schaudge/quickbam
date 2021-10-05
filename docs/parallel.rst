Parallel processing of BAM files
================================

The main advantage of using libmmbam is to achieve stunningly fast speed via
parallel data processing, with the help of the BAM index file. Many data
processing tasks can be expressed in a map-reduce fashion where "map"
operations can be applied independently to many genomic regions, followed by a
"reduce" operation to combine the partial results into a final result.

In this part of the documentation, we will illustrate the basic constructs for
implementing such a parallelization pattern.

Parallelize by regions
======================

If a computation task is based on reads, it is often parallelizable over many
regions. Such an example (albeit trivial) is to count the number of reads in a
BAM file. A program can go through all the reads one by one from the beginning
to the end, or, it can process many small, non-overlapping regions in parallel,
and then sum the counts of each region up. Here we demonstrate how to take the
second approach, and parallelize on each 16kb genomic window, which maps
directly to the BAM linear index. Libmmbam offers a convinence function, 
index2region, to return a vector of (start, end] virtual offset pairs, which
can be used to load the BAM reads of the corresponding regions. We use OpenMP
as the parallelization mechanism for this example

.. code-block:: cpp

   #include <iostream>
   #include <fstream>
   #include <mfile.h>
   #include <bam.h>
   #include <index.h>

   int main(int argc, const char *const argv[]) {

       // we ignore bound and error checking in the example

       auto bam_filename = argv[1];
       auto bai_filename = std::string(argv[1] + ".bai");

       auto mfile = mfile_open(bam_filename);
       auto index = index_read(std::ifstream(bai_filename));

       uint64_t read_count = 0;

       // generate parallelization regions
       auto regions = index_to_regions(index, mfile->size);

       #pragma omp parallel for reduction(+:read_count)
       for(size_t i=0; i<regions.size(); i++) {
           auto bam_records = bam_load_block(mfile, regions[i].first, regions[i].second);
           auto block_count = bam_count_records(bam_records);
           read_count += block_count;
       }

       std::cout<<"total records="<<read_count<<std::endl;
       index_free(index);
       return 0;
   }

Parallelize by positions (mpileup)
==================================

Another pattern of compute tasks is to gather all reads overlapping a specific
position, over many independent positions. If the task is to process each read
independently, it can be implemented in a similar fashion as the example above.
However, if it is the bases from all these reads overlapping the position that
need to be processed, a pileup engine is then necessary. Libmmbam provides a
multiple input pileup engine, which can be invoked in parallel over many
positions. Here we provide some simple concepts on how to use the pileup
engine. For a full working program example, please see the snp-pileup-tbb
program in the "code_example" directory of the repository.

.. code-block:: cpp

   #include <mfile.h>
   #include <bam.h>
   #include <index.h>
   #include <mpileup.h>

   #include <fstream>

   #include <tbb/parallel_for_each.h>
   
   struct work_items {
       // location to pileup
       int32_t ref_id, pos;
       // the reference base at this location
       uint8_t ref_base;

       // results after processing the piled up location
       // This is necessary because parallel jobs will need independent
       // locations to write their results, else a mutex will be necessary
       int ref_count[2];    // count ref alleles per input buffer
       int other_count[2];  // count other alleles per input buffer
   }

   int main(void) {
       // mpileup engine supports multiple input files.
       // This example shows how to pileup two files simultaneously

       auto mfile1 = mfile_open("bam_file1.bam");
       auto mfile2 = mfile_open("bam_file2.bam");

       auto index1 = index_read(std::ifstream("bam_file1.bam.bai"));
       auto index2 = index_read(std::ifstream("bam_file2.bam.bai"));

       mfiles_t mfiles{mfile1, mfile2}
       indices_t indices{index1, index2}

       // the following line calls some function to generate parallel work items
       // which should fill the ref_id, pos, and ref_base fields of each item,
       // and initialize the ref_count and other_count fields to 0
       std::vector<work_items> parallel_items = work_item_gen_func();

       tbb::parallel_for_each(
           parallel_items.cbegin(),
           parallel_items.cend(),
           [&](auto& r) {

               // the filter lambda returns false if the mapping quality of
               // a read is below a given a hardcoded threshold (1)
               auto filter_predicate = [](const auto& bam_rec) {
                   if(bam_rec.mapq<1) return false;
                   return true;
               }

               // The visitor function lambda will be called at each piled up
               // location with a mpileup_t struct as the parameter
               auto visitor_func = [&r](const auto& p) {

                   // skip all positions before the desired position
                   if(p.pos < r.pos) return true;

                   // halt pileup if we are past the desired position
                   if(p.pos > r.pos) return false;

                   // here implies that the piled up position is equal to r.pos
                   // count reference and other reads per input buffer
                   for(size_t i_file = 0; i_file < 2; i_file++) {
                       auto* buffer = p.reads_buffer[i_file]->data();

                       // iterate over reads in buffer i_file
                       for(size_t i=0; i<p.depth[i_file]; i++) {
                           auto& info = p.get_info(i_file, i);
                           const bam_rec_t* bam_record = BAMREF(buffer + info.offset);

                           // only count non-deletion alleles
                           if(!info.is_deletion) {
                               if(bam_bqual_ptr(bam_record)[info.qpos] < 1)
                                   continue; // ignore low quality bases
                               auto seq = bam_seq_ptr(bam_record);
                               auto base = bam_unpack_base(seq, info.qpos);
                               if(base == r.ref_base) r.ref_count[i_file]++;
                               else r.other_count[i_file]++;
                           }
                       } // end for-each read
                   } // end for-each buffer

                   return true;
               }

               // call mpileup engine
               mpileup(mfiles, indices, r.ref_id, r.pos, r.pos+1,
                       filter_predicate, visitor_func);

       });

       // parallel_items now contain the results
       // post-processing
       // output

       index_free(index1);
       index_free(index2);

       return 0;
   }

