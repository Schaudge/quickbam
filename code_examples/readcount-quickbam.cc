#include<iostream>
#include<fstream>

#include "quickbam/mfile.h"
#include "quickbam/bam.h"
#include "quickbam/index.h"
#include "quickbam/slicer.h"

// perform parallel read of a bam file.

int main(int argc, const char *const argv[]) {

    if(argc < 2) {
        std::cout<<"usage: "<<argv[0]<<" <bam-file> [bai-file]"<<std::endl;
        exit(0);
    }

    auto bam_filename = argv[1];
    auto bai_filename = argc < 3 ? std::string(argv[1]) + ".bai" : std::string(argv[2]);

    file_slicer_t data(bam_filename);
    auto index = index_read(std::ifstream(bai_filename));

    uint64_t records = 0;

    // create iteration intervals
    auto regions = index_to_regions(index, data.size());


#pragma omp parallel for reduction(+:records)
    for(size_t i=0; i<regions.size(); i++) {
        if(regions[i].first == regions[i].second) throw std::runtime_error("same region start and end");
        auto bam_records = bam_load_block(data, regions[i].first, regions[i].second);
        auto record_count = bam_count_records(bam_records);
        records += record_count;
    }

    std::cout<<"total records="<<records<<std::endl;

    index_free(index);
    return 0;
}
