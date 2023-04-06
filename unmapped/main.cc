#include <stdio.h>
#include <fstream>
#include <cmath>
#include <omp.h>

#include <quickbam/bam.h>
#include <quickbam/slicer.h>

// 1. Start with slicer
// 2. Divide it into N chunks based on parallelism
// 3. 

using namespace std;

#define PRINT(x) std::cout << #x << ": " << x << std::endl


void print_rec(const bam_rec_t& rec) {
    cout
        << "block_size: " << rec.block_size
        << endl
        << "ref_id: " << rec.ref_id
        << endl
        << "pos: " << rec.pos
        << endl
        << "l_read_name: " << (int)rec.l_read_name
        << endl
        << "l_seq: " << rec.l_seq
        << endl << endl;
}

void print_region(const region& region) {
    cout
        << index_coffset(region.first)
        << "-"
        << index_coffset(region.second)
        << endl;
}

void print_block(const bgzf_block_t* b) {
    cout
        << "id1: " << (int)b->id1 << endl
        << "id2: " << (int)b->id2 << endl
        << "si1: " << (int)b->si1 << endl
        << "si2: " << (int)b->si2 << endl
        << "bsize: " << b->bsize  << endl
        << "isize: " << bgzf_isize(*b) << endl
        << endl;
}











int main(int argc, char** argv) {

    file_slicer_t slicer(argv[1]);

    auto regions = bam_to_regions(slicer);

    uint64_t records = 0;

#pragma omp parallel for reduction(+:records)
    for(size_t i=0; i<regions.size(); i++) {
        if(regions[i].first == regions[i].second) throw std::runtime_error("same region start and end");
        auto bam_records = bam_load_block(slicer, regions[i].first, regions[i].second);
        //cout << regions[i].first << endl;
        auto record_count = bam_count_records(bam_records);
        records += record_count;
    }

    std::cout<<records<<std::endl;

    return 0;
}
