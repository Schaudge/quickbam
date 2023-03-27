#include "quickbam/mbgzf.h"
#include "quickbam/bam.h"
#include "quickbam/index.h"
#include <vector>
#include <iostream>
#include <cassert>

#include <tbb/parallel_for.h>


bool bam_buffer_contains_header(const std::vector<uint8_t>& buffer) {
    // fixed part of the header
    size_t curr_size = 8;
    if(buffer.size() < curr_size) return false;
    const uint8_t *buffptr = &buffer[0];

    // text and n_ref
    uint32_t l_text = *(uint32_t*)(buffptr+4);
    curr_size += (l_text + 4);
    if(buffer.size() < curr_size) return false;

    uint32_t n_ref = *(uint32_t*)(buffptr + l_text + 8);
    for(uint32_t i_ref = 0; i_ref<n_ref; ++i_ref) {
        // need l_name
        curr_size += 4;
        if(buffer.size() < curr_size) return false;
        uint32_t l_name = *(uint32_t *)(buffptr + curr_size - 4);
        // need name and l_ref
        curr_size += l_name + 4;
        if(buffer.size() < curr_size) return false;
    }
    return true;
};

size_t bam_count_records(const std::vector<uint8_t>& buffer) {
    auto bam_it = bam_iterator(reinterpret_cast<const bam_rec_t*>(&buffer[0]));
    auto bam_it_beg = bam_iterator(reinterpret_cast<const bam_rec_t*>(&buffer[0]));
    auto bam_it_end = bam_iterator(reinterpret_cast<const bam_rec_t*>(&buffer[0] + buffer.size()));
    size_t total_recs = 0;
    while(bam_it < bam_it_end) {
        total_recs++;
        bam_it++;
    }

    return total_recs;
}

int16_t bam_query_length(const bam_rec_t* b) {
    uint32_t ql = 0;
    const uint32_t *cigar_ops = (const uint32_t *)(
            (const uint8_t *)(b)
            + sizeof(bam_rec_t)
            + sizeof(char) * b->l_read_name );
    for(int op = 0; op < b->n_cigar_op; op++) {
        switch(cigar_ops[op] & 0xF) {
            case 0: // M
            case 2: // D
            case 7: // =
            case 8: // X
                ql += cigar_ops[op] >> 4;
                break;
            default:
                break;
        }
    }
    return ql;
}
