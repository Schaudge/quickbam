#include "quickbam/mbgzf.h"
#include "quickbam/bam.h"
#include "quickbam/index.h"
#include <vector>
#include <iostream>
#include <cassert>

#include <tbb/parallel_for.h>

#include <taskflow/taskflow.hpp>



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

struct bam_header_aux_data_t {
    char tag[2];
    char val_type;
    uint8_t vardata[];
};

int32_t bam_consume_aux_item(const uint8_t* ptr, size_t max_size) {
    auto aux = reinterpret_cast<const bam_header_aux_data_t*>(ptr);

    int32_t consumed = sizeof(bam_header_aux_data_t);

    switch (aux->val_type) {
        case 'A':
        case 'c':
        case 'C': {
            consumed += 1;
            break;
        }
        case 's':
        case 'S': {
            consumed += 2;
            break;
        }
        case 'i':
        case 'I':
        case 'f': {
            consumed += 4;
            break;
        }
        case 'H':
        case 'Z': {
            bool found = false;
            for (size_t i = 0; i < max_size; i++) {
                consumed += 1;
                if (aux->vardata[i] == '\0') {
                    found = true;
                    break;
                }
            }

            if (!found) {
                return -1;
            }
            break;
        }
        case 'B': {
            char subtype = aux->vardata[0];
            consumed += 1;

            uint32_t count = *reinterpret_cast<const uint32_t*>(&aux->vardata[1]);
            consumed += 4;

            switch(subtype) {
                case 'c':
                case 'C': {
                    consumed += count;
                    break;
                }
                case 's':
                case 'S': {
                    consumed += (2 * count);
                    break;
                }
                case 'i':
                case 'I':
                case 'f': {
                    consumed += (4 * count);
                    break;
                }
            }
            break;
        }
        default: {
            return -1;
            break;
        }
    }

    return consumed;
}


static bool bam_valid_aux_data(const bam_rec_t* r) {
    uint64_t aux_offset =
        //sizeof(bam_rec_t) +
        r->l_read_name * sizeof(char) + // read name
        r->n_cigar_op * sizeof(uint32_t) + // cigar
        ((r->l_seq + 1) / 2) * sizeof(uint8_t) + // seq
        r->l_seq * sizeof(char) // qual
    ;

    auto aux_ptr = &r->vardata[aux_offset];

    size_t total_consumed = sizeof(bam_rec_t) - sizeof(uint32_t) + aux_offset;

    for (;;) {
        
        auto consumed = bam_consume_aux_item(aux_ptr, r->block_size);
        if (consumed == -1) {
            return false;
        }

        aux_ptr += consumed;
        total_consumed += consumed;

        if (total_consumed == r->block_size) {
            break;
        }
        else if (total_consumed > r->block_size) {
            return false;
        }
    }

    return true;
}

bool bam_is_valid(bam::Header& header, const bam_rec_t* r) {

    auto& refs = header.refs();
    auto n_ref = static_cast<int32_t>(refs.size());

    if (!(-1 <= r->ref_id && r->ref_id < n_ref)) {
        return false;
    }
    if (!(-1 <= r->next_ref_id && r->next_ref_id < n_ref)) {
        return false;
    }

    if (r->ref_id == -1) {
        if (r->pos != 0) {
            //return false;
        }
    }
    else {
        if (!(r->pos < refs[r->ref_id].l_ref)) {
            return false;
        }
    }

    if (r->next_ref_id == -1) {
        if (r->next_pos != 0) {
            //return false;
        }
    }
    else {
        if (!(r->next_pos < refs[r->next_ref_id].l_ref)) return false;
    }

    auto read_name = r->vardata;

    if (read_name[r->l_read_name - 1] != '\0') {
        return false;
    }


    for (size_t i = 0; i < r->l_read_name - 1; i++) {
        if (read_name[i] < '!' || read_name[i] == '@' || read_name[i] > '~') {
            return false;
        }
    }

    //if (!(r->l_seq < r->block_size)) return false;

    auto valid_aux_data = bam_valid_aux_data(r);
    if (!valid_aux_data) return false;

    return true;
}
