#include <stdio.h>
#include <fstream>
#include <cmath>

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

bool is_bam(const bam_rec_t* r) {
    return r->block_size > 100 &&
        r->block_size < 10000 &&
        r->l_seq > 1 &&
        r->l_seq < 500;
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
std::vector<region> slicer_to_regions(SLICER_T slicer, size_t start_offset, size_t end_offset) {

    std::vector<region> regions;

    const size_t CHUNK_SZ = 1024*1024;

    uint64_t region_start = 0;
    uint64_t ioffset = 0;
    bool first = true;
    bool found = false;

    //for (size_t start_idx = 0; start_idx < 10000000; start_idx += CHUNK_SZ) {
    for (size_t start_idx = start_offset; start_idx < end_offset; start_idx += CHUNK_SZ) {
        //cout << "start_idx: " << start_idx << endl;
        auto block_idx = bgzf_find_next_block(slicer, start_idx);

        bgzf_slicer_iterator_t bgzf_it(slicer, block_idx);
        auto bgzf_end = bgzf_it.end();

        uint64_t coffset = block_idx;

        found = false;
        while (bgzf_it != bgzf_end) {

            const bgzf_block_t& block = *bgzf_it;

            auto block_bytes = bgzf_inflate(block);

            for (size_t uoffset = 0; uoffset < block_bytes.size(); uoffset++) {
                // TODO: apparently with HG002.GRCh38.2x250.bam j is never greater
                // than 0?
                auto bam_rec = reinterpret_cast<const bam_rec_t*>(&block_bytes[uoffset]);

                if (is_bam(bam_rec)) {
                    found = true;
                    ioffset = calc_ioffset(coffset, uoffset);
                    break;
                }
            }

            if (found) {
                break;
            }

            coffset += (block.bsize + 1);
            bgzf_it++;
        }

        if (found) {
            if (first) {
                first = false;
                region_start = ioffset;
            }
            else {
                regions.push_back({ region_start, ioffset });
                region_start = ioffset;
            }
        }

        assert(found);
    }


    if (index_coffset(region_start) < end_offset) {
        regions.push_back({ region_start, calc_ioffset(end_offset, 0) });
    }

    return regions;
}


template<typename SLICER_T>
std::vector<region> slicer_to_regions(SLICER_T slicer, size_t start_offset) {
    return slicer_to_regions(slicer, start_offset, slicer.size());
}

template<typename SLICER_T>
std::vector<region> slicer_to_regions(SLICER_T slicer) {
    return slicer_to_regions(slicer, 0, slicer.size());
}


int main(int argc, char** argv) {

    file_slicer_t slicer(argv[1]);

    //uint64_t idx = stoull(argv[2]);
    //cout << "idx:" << idx << endl;
    //auto slice = slicer.slice(idx, idx + BGZF_MAX_BLOCK_SIZE);
    //auto blk = reinterpret_cast<const bgzf_block_t*>(slice.get());

    //print_block(blk);

    //auto buf = bgzf_inflate(*blk);

    //auto bam_read = reinterpret_cast<const bam_rec_t*>(buf.data());

    //print_rec(*bam_read);

    //auto regions = slicer_to_regions(slicer);



    auto index = index_read(std::ifstream(std::string(argv[1]) + ".bai"));
    auto index_regions = index_to_regions(index, slicer.size());

    auto last_mapped_region = index_regions[index_regions.size() - 1];


    cout << "slicer_to_regions" << endl;
    auto regions = slicer_to_regions(slicer);
    //auto regions = slicer_to_regions(slicer, 0, 4733653485);
    //auto regions = slicer_to_regions(slicer, 1000000);
    
    for (auto& region : regions) {
        print_region(region);
    }

    uint64_t min_size = 1000000000;
    uint64_t max_size = 0;
    vector<uint64_t> sizes;
    for (auto& region : regions) {

        //print_region(region);

        auto size = index_coffset(region.second) - index_coffset(region.first);
        sizes.push_back(size);
        if (size < min_size) {
            min_size = size;
        }
        if (size > max_size) {
            max_size = size;
        }

        //if (size > 5000000) {
        //    cout << size << endl;
        //}
    }

    if (sizes.size() > 0) {
        uint64_t sum = 0;
        for (auto& size : sizes) {
            sum += size;
        }
        auto avg = sum / sizes.size();

        cout << "number of regions: " << regions.size() << endl;
        cout << "min region size: " << min_size << endl;
        cout << "max region size: " << max_size << endl;
        cout << "avg region size: " << avg << endl;
    }

    uint64_t records = 0;

    cout << "counts" << endl;

#pragma omp parallel for reduction(+:records)
    for(size_t i=0; i<regions.size(); i++) {
        if(regions[i].first == regions[i].second) throw std::runtime_error("same region start and end");
        auto bam_records = bam_load_block(slicer, regions[i].first, regions[i].second);
        //cout << regions[i].first << endl;
        auto record_count = bam_count_records(bam_records);
        records += record_count;
    }

    std::cout<<"total records="<<records<<std::endl;



    




    index_free(index);

    return 0;
}





//struct byte_result_t {
//    uint8_t value;
//    bool error;
//};
//
//const size_t BYTE_SLICER_CHUNK_SIZE = 64*1024;
//
//template<typename SLICER_T>
//class byte_slicer_t {
//    private:
//        SLICER_T slicer;
//        typename SLICER_T::ptr_t cur_slice;
//        size_t cur_offset;
//
//        void update_slice() {
//            cur_slice = slicer.slice(cur_offset, cur_offset + BYTE_SLICER_CHUNK_SIZE);
//        }
//
//    public:
//        byte_slicer_t(SLICER_T slicer) : slicer(slicer), cur_offset(0) {
//            update_slice();
//        }
//
//        byte_result_t get_at(size_t index) {
//            byte_result_t result;
//            
//            if (index >= slicer.size()) {
//                result.error = true;
//                return result;
//            }
//            else {
//                result.error = false;
//            }
//
//            if (index < cur_offset || index >= (cur_offset + BYTE_SLICER_CHUNK_SIZE)) {
//                cur_offset = index;
//                update_slice();
//            }
//
//            auto adjusted_index = index - cur_offset;
//            result.value = cur_slice[adjusted_index];
//
//            return result;
//        }
//
//        size_t size() {
//            return slicer.size();
//        }
//};

//ofstream out_file("data.tsv");
//
//    out_file << "block_size\tref_id\tpos\tl_read_name\tl_seq" << endl;
//
//    size_t record_counter = 0;
//    size_t diff_counter = 0;
//    uint32_t l_seq = 0;
//    auto first = true;
//
//    for (size_t i = 1; i < regions.size(); i++) {
//    //for (size_t i = 1; i < 2; i++) {
//
//        auto region = regions[i];
//
//        auto bam_records = bam_load_block(slicer, region.first, region.second); 
//        auto bam_it_beg = bam_iterator(reinterpret_cast<const bam_rec_t*>(&bam_records[0]));
//        auto bam_it_end = bam_iterator(reinterpret_cast<const bam_rec_t*>(&bam_records[0] + bam_records.size()));
//        auto bam_it = bam_it_beg;
//        while(bam_it < bam_it_end) {
//            auto r = *bam_it;
//
//            record_counter++;
//
//            if (first) {
//                l_seq = r.l_seq;
//                first = false;
//            }
//
//            uint32_t diff = 0;
//            if (r.l_seq > l_seq) {
//                diff = r.l_seq - l_seq;
//            }
//            else {
//                diff = l_seq - r.l_seq;
//            }
//
//            if (diff > 1) {
//                diff_counter++;
//                cout << diff << endl;
//                printRec(r);
//            }
//            //out_file
//            //    << r.block_size << "\t"
//            //    << r.ref_id << "\t"
//            //    << r.pos << "\t"
//            //    << (int)r.l_read_name << "\t"
//            //    << r.l_seq << endl;
//            bam_it++;
//        }
//    }
//
//    cout << "record_counter: " << record_counter << endl;
//    cout << "diff_counter: " << diff_counter << endl;
//
//
//    //auto bam_block = bam_load_block(slicer, regions[1].first, regions[1].second); 
//    //auto rec = BAMREF(bam_block.data());
//    //printRec(*rec);
//
//    //for (size_t i = 0; i < 3; i++) {
//    //    rec = BAM_NEXT(rec);
//    //    printRec(*rec);
//    //}
//
//    out_file.close();
//
//
