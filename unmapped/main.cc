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

bool bam_is_valid(const bam_rec_t* r) {
    bool basic_checks = r->l_seq > 0 &&
        r->l_seq < 500 &&
        r->block_size > r->l_seq &&
        r->block_size < 2000 &&
        -1 <= r->ref_id &&
        -1 <= r->next_ref_id &&
        //r->next_ref_id < r->n_ref &&
        r->l_read_name + 4*r->n_cigar_op + r->l_seq < r->block_size;

    if (!basic_checks) return false;

    if (bam_read_name(r)[r->l_read_name - 1] != '\0') return false;

    return true;
}


const uint32_t BGZF_MAX_BLOCK_SIZE = 64*1024;

template<typename SLICER_T>
size_t bgzf_find_next_block(SLICER_T slicer, size_t start_offset) {

    // A valid BGZF file should never have more than 64k bytes between blocks.
    // Therefore when searching for a valid block header (bgzf_block_t), the
    // worst case is that the block starts at 64k-1 from the start of the
    // search. So we should never need more than 64k+sizeof(bgzf_block_t)-1
    // bytes to find a valid header.
    const uint32_t SLICE_SIZE = BGZF_MAX_BLOCK_SIZE + sizeof(bgzf_block_t) - 1;

    auto slice = slicer.slice(start_offset, start_offset + SLICE_SIZE);

    for (size_t i = 0; i < SLICE_SIZE; i++) {
        auto blk = reinterpret_cast<const bgzf_block_t*>(&(slice.get()[i]));
        if (blk->id1 == 31 && blk->id2 == 139 && blk->si1 == 66 && blk->si2 == 67 && blk->slen == 2) {
            return start_offset + i;
        }
    }

    // TODO: properly handle error conditions
    assert(0 != 1);
    return 0;
}

uint64_t calc_ioffset(uint64_t coffset, uint64_t uoffset) {
    return (coffset << 16) | uoffset;
}


template<typename SLICER_T>
std::vector<region> bam_to_regions(SLICER_T slicer, size_t start_offset, size_t end_offset) {

    const size_t CHUNK_SZ = 1024*1024;

    auto num_chunks = (end_offset - start_offset) / CHUNK_SZ;
    
    std::vector<size_t> region_starts(num_chunks);

#pragma omp parallel for
    for (size_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {

        auto start_idx = chunk_idx*CHUNK_SZ;

        auto block_idx = bgzf_find_next_block(slicer, start_idx);

        bgzf_slicer_iterator_t bgzf_it(slicer, block_idx);
        auto bgzf_end = bgzf_it.end();

        uint64_t coffset = block_idx;

        bool found = false;
        while (bgzf_it != bgzf_end) {

            const bgzf_block_t& block = *bgzf_it;

            auto block_bytes = bgzf_inflate(block);

            for (size_t uoffset = 0; uoffset < block_bytes.size(); uoffset++) {
                // TODO: apparently with HG002.GRCh38.2x250.bam j is never greater
                // than 0?
                auto bam_rec = reinterpret_cast<const bam_rec_t*>(&block_bytes[uoffset]);

                if (bam_is_valid(bam_rec)) {
                    found = true;
                    uint64_t ioffset = calc_ioffset(coffset, uoffset);
                    region_starts[chunk_idx] = ioffset;
                    break;
                }
            }

            if (found) {
                break;
            }

            coffset += (block.bsize + 1);
            bgzf_it++;
        }

        assert(found);
    }

    std::vector<region> regions;

    auto first = true;
    uint64_t region_start = 0;
    for (auto& cur_start : region_starts) {
        if (first) {
            first = false;
        }
        else {
            auto region_end = cur_start;
            regions.push_back({ region_start, region_end });
        }

        region_start = cur_start;
    }

    auto last_start = region_starts[num_chunks - 1];

    if (index_coffset(last_start) < end_offset) {
        regions.push_back({ last_start, calc_ioffset(end_offset, 0) });
    }

    return regions;
}


template<typename SLICER_T>
std::vector<region> bam_to_regions(SLICER_T slicer, size_t start_offset) {
    return bam_to_regions(slicer, start_offset, slicer.size());
}

template<typename SLICER_T>
std::vector<region> bam_to_regions(SLICER_T slicer) {
    return bam_to_regions(slicer, 0, slicer.size());
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
