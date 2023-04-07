#ifndef BAM_H
#define BAM_H

/****************************************************************
 * BAM file structure APIs implemented on memory buffers
*****************************************************************/

#include <stdint.h>
#include <vector>
#include "nfo_iterator.h"
#include "mfile.h"
#include "index.h"
#include "quickbam/mbgzf.h"
#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>

#define BYTEREF(b) ((const uint8_t*)(b))
#define BAMREF(b)  ((const bam_rec_t *)(b))
#define BAM_NEXT(b) (BAMREF(BYTEREF(b) + (b)->block_size + 4))
#define READ_IN_REGION(b, rbegin, rend) ((b)->pos < (rend) && (b)->pos + bam_query_length(b) >= (rbegin))

#define READ_PAIRED 0x1
#define READ_MAPPED_PROPER_PAIR 0x2
#define READ_UNMAPPED 0x4
#define MATE_UNMAPPED 0x8
#define READ_REVERSE_STRAND 0x10
#define MATE_REVERSE_STRAND 0x20
#define READ_FIRST_IN_PAIR 0x40
#define READ_SECOND_IN_PAIR 0x80
#define READ_SECONDARY_ALIGNMENT 0x100
#define READ_FAILED_QUALITY_CHECKS 0x200
#define READ_PCR_OPTICAL_DUPLICATE 0x400
#define READ_SUPPLEMENTARY 0x800
#define READ_PRIMARY 0x900

#define CHUNK_SIZE 16384

//! This struct mirrors the byte layout of uncompressed bam header
struct bam_header_t {
    char magic[4];      // BAM\1
    uint32_t l_text;    // Length of header text
    uint8_t vardata[];  // Data of variable length

    // char text[l_text];   // Header text
    // uint32_t n_ref;      // Number of reference sequences
    // struct[] {
    //   uint32_t l_name;   // Length of the reference name
    //   char name[l_name]; // Reference name
    //   uint32_t l_ref;    // Length of the reference sequence
    // }
};

//! this struct mirrors the byte layout of uncompressed bam records
struct bam_rec_t {
    uint32_t block_size;    /**< Total length of alignment record except this field */
    int32_t  ref_id;        /**< Reference sequence id; -1 if unmapped */
    int32_t  pos;           /**< 0-based left-most position */
    uint8_t  l_read_name;   /**< Length of read_name */
    uint8_t  mapq;          /**< Mapping quality */
    uint16_t bin;           /**< BAI index bin */
    uint16_t n_cigar_op;    /**< Number of operations in CIGAR */
    uint16_t flag;          /**< Bitwise flags */
    uint32_t l_seq;         /**< Length of SEQ */
    int32_t  next_ref_id;   /**< Reference id of the next segment */
    int32_t  next_pos;      /**< 0-based left-most position of the next segment */
    int32_t  tlen;          /**< Template length */
    uint8_t  vardata[];     /**< Data of variable length */
    
    //  char[l_read_name]       // Read name
    //  uint32_t[n_cigar_op]    // CIGAR
    //  uint8_t[(l_seq + 1)/2]  // 4-bit encoded read seq
    //  char[l_seq]             // Phred-scale base quality

    //  struct[] {              // TAGS
    //      char[2] tag;        // Two-character tag
    //      char val_type;      // Value type
    //      value;              // Value (type by val_type)
};

using bam_iterator = nfo_iterator<bam_rec_t, uint32_t, 0, sizeof(uint32_t)>;

//! check if the given buffer contains a complete bam header record
bool bam_buffer_contains_header(const std::vector<uint8_t>& buffer); 

//! returns the number of complete bam records found in the given buffer
size_t bam_count_records(const std::vector<uint8_t>& buffer);

//! returns the total length of the query string of the bam record
int16_t bam_query_length(const bam_rec_t* b);

bool bam_is_valid(const bam_rec_t* r);

//! returns the read name of the bam record
inline std::string bam_read_name(const bam_rec_t* b) {
    return std::string(b->vardata, b->vardata + b->l_read_name - 1);
}

//! returns the pointer to the cigar operations array of the bam record
inline const uint32_t* bam_cigar_ptr(const bam_rec_t* b) {
    return reinterpret_cast<const uint32_t *>(b->vardata
        + sizeof(char) * b->l_read_name);
}

//! returns the pointer to the query sequence string of the bam record
inline const uint8_t* bam_seq_ptr(const bam_rec_t* b) {
    return b->vardata
        + sizeof(char) * b->l_read_name
        + sizeof(uint32_t) * b->n_cigar_op;
}

//! returns the pointer to the base quality array of the bam record
inline const char* bam_bqual_ptr(const bam_rec_t* b) {
    return reinterpret_cast<const char *>(b->vardata
            + sizeof(char) * b->l_read_name
            + sizeof(uint32_t) * b->n_cigar_op
            + sizeof(uint8_t) * ((b->l_seq + 1) / 2));
}

//! extract the sequence from a 4-bit compact representation sequence array at given position
inline uint8_t bam_unpack_base(const uint8_t* seq, uint32_t qpos) {
    if(qpos % 2 == 0) // even base, higher 4 bits
        return seq[qpos / 2] >> 4;
    else // odd base, lower 4 bits
        return seq[qpos / 2] & 0xF;
}

const uint8_t bam_a_lo = 0x01;
const uint8_t bam_c_lo = 0x02;
const uint8_t bam_g_lo = 0x04;
const uint8_t bam_t_lo = 0x08;
const uint8_t bam_a_hi = 0x10;
const uint8_t bam_c_hi = 0x20;
const uint8_t bam_g_hi = 0x40;
const uint8_t bam_t_hi = 0x80;


inline char byte2base_lo(const uint8_t base) {
    return base & bam_a_lo ? 'A' :
        base & bam_c_lo ? 'C' :
        base & bam_g_lo ? 'G' :
        base & bam_t_lo ? 'T' : 'N';
};

inline char byte2base_hi(const uint8_t base){
    return base & bam_a_hi ? 'A' :
        base & bam_c_hi ? 'C' :
        base & bam_g_hi ? 'G' :
        base & bam_t_hi ? 'T' : 'N';
};

//! loads a block of bam records from a region
//! \param data Slicer object which abstracts over the source BAM file.
//! \param ioffset_first the ioffset of the beginning of the region. ioffset is the virtual offset defined per SAM file specification, section 4.1.1
//! \param ioffset_last the ioffset of the end of the region.
//! \return a byte vector of all bam files contained in the region
template<class SLICER_T>
std::vector<uint8_t> bam_load_block(SLICER_T data, uint64_t ioffset_first, uint64_t ioffset_last) {

    if(ioffset_first == ioffset_last) return std::vector<uint8_t>();

    auto coffset_first = index_coffset(ioffset_first);
    auto uoffset_first = index_uoffset(ioffset_first);

    auto coffset_last = index_coffset(ioffset_last);
    auto uoffset_last = index_uoffset(ioffset_last);

    // TODO: We're adding BGZF_MAX_BLOCK_SIZE here because coffset_last
    // represents the *beginning* of the last block. BGZF_MAX_BLOCK_SIZE
    // guarantees we slice enough data to decompress the whole block. It might
    // be better to read the block header and figure out exactly how much we
    // need to read, but I'm guessing the extra IO might actually be slower.
    auto load_end = coffset_last + BGZF_MAX_BLOCK_SIZE;
    auto slice = data.slice(coffset_first, load_end < data.size() ? load_end : data.size());

    auto bgzf_block_first = reinterpret_cast<const bgzf_block_t*>(slice.get());
    auto bgzf_block_last = reinterpret_cast<const bgzf_block_t*>(slice.get() + (coffset_last - coffset_first));

    // special case: coffsets are the same
    if(coffset_first == coffset_last) {
        auto inflated_bytes = bgzf_inflate(*bgzf_block_first);

        auto first_byte = inflated_bytes.cbegin() + uoffset_first;
        auto last_byte  = uoffset_last == uoffset_first ? 
            inflated_bytes.cend() :
            inflated_bytes.cbegin() + uoffset_last;

        return std::vector<uint8_t>(first_byte, last_byte);
    }


    // inflate from bgzf_block_first to bgzf_block_last
    // trim uoffset_first from the beginning
    // trim uoffset_last from the end

    //auto bgzf_block_it = bgzf_block_first;

    /**** SEQUENCIAL DECOMPRESSION ****/
    //std::vector<uint8_t> inflated_bytes_s;
    //inflated_bytes_s = bgzf_inflate_range(
    //        begin<const uint8_t>(mfile) + coffset_first,
    //        coffset_last - coffset_first);

    /**** PARALLEL DECOMPRESSION ****/
    auto inflated_bytes = bgzf_inflate_range_p(
            slice.get(),
            coffset_last - coffset_first, [](
        auto* src, auto src_len, auto dest_len,
        auto& src_off_vector, auto& dest_off_vector,
        auto inflate) -> std::vector<uint8_t> {

        std::vector<uint8_t> buffer;
        buffer.resize(dest_len);

        tbb::this_task_arena::isolate([&] {

            tbb::parallel_for(tbb::blocked_range<size_t>(0, src_off_vector.size()), [&](const auto& r){

                auto sl = r.end() == src_off_vector.size() ? 
                src_len                 - src_off_vector[r.begin()] : 
                src_off_vector[r.end()] - src_off_vector[r.begin()];

                auto dl = r.end() == dest_off_vector.size() ? 
                dest_len                 - dest_off_vector[r.begin()] :
                dest_off_vector[r.end()] - dest_off_vector[r.begin()];

                inflate(src+src_off_vector[r.begin()], sl, &buffer[dest_off_vector[r.begin()]], dl);

            });
        });

        

        return buffer;
    });

    if(coffset_last < data.size()) {
        auto last_block = bgzf_inflate(*bgzf_block_last);
        inflated_bytes.insert(inflated_bytes.end(), last_block.cbegin(), last_block.cbegin() + uoffset_last);
    }

    return std::vector<uint8_t>(inflated_bytes.cbegin() + uoffset_first, inflated_bytes.cend());
}

template<class T, class V>
bool extend_buffer(T& bgzf_it, T& bgzf_end, V& buffer, size_t min_size) {
    if(buffer.size() >= min_size) return true;
    if(bgzf_it == bgzf_end) return false;
    while(bgzf_it != bgzf_end && buffer.size() < min_size) {
        auto inflated = bgzf_inflate(*bgzf_it);
        bgzf_it++;
        buffer.insert(buffer.end(), inflated.cbegin(), inflated.cend());
    }
    return buffer.size() >= min_size;
}


//! Helper function to load all bam records from a specific genomic region
//
//! \param data Slicer object which abstracts over the source BAM file.
//! \param index The index created on the BAM file
//! \param reg_id The numerical id of the reference contig
//! \param region_start The genomic coordinate of the start of the region
//! \param region_end The genomic coordinate of the end of the region
//! \return A byte vector containing the bam records within the specified region
template<class SLICER_T>
std::vector<uint8_t> bam_load_region(SLICER_T data, const index_t& index, int32_t ref_id, int32_t region_start, int32_t region_end) {

    if(ref_id >= index.n_ref) throw std::runtime_error("reference not found");
    if(index.ref[ref_id].n_intv == 0) throw std::runtime_error("empty reference");

    using intv_idx_t = decltype(index.ref[ref_id].n_intv);

    // start_interval
    intv_idx_t start_intv = region_start / CHUNK_SIZE;
    intv_idx_t end_intv   = region_end   / CHUNK_SIZE;

    if(start_intv > index.ref[ref_id].n_intv - 1)
        start_intv = index.ref[ref_id].n_intv - 1;

    while(end_intv < index.ref[ref_id].n_intv && 
            index.ref[ref_id].ioffset[start_intv] == index.ref[ref_id].ioffset[end_intv])
        end_intv++;

    if(end_intv > index.ref[ref_id].n_intv - 1)
        end_intv = index.ref[ref_id].n_intv - 1;

    // load bgzf block in batch
    std::vector<uint8_t> bam_buffer_preload = bam_load_block(
            data,
            index.ref[ref_id].ioffset[start_intv],
            index.ref[ref_id].ioffset[end_intv]);

    auto bam_iter = reinterpret_cast<const bam_rec_t *>(bam_buffer_preload.data());
    auto buffer_end = bam_buffer_preload.data() + bam_buffer_preload.size();
    size_t first_in_region = -1;

    if(bam_buffer_preload.size() > 0) {
        auto bam_next = BYTEREF(bam_iter) + bam_iter->block_size + 4;
        // reads returned. See if the region is contained
        while(bam_next < buffer_end) {
            if( READ_IN_REGION(bam_iter, region_start, region_end) && first_in_region == -1) {
                first_in_region = BYTEREF(bam_iter) - bam_buffer_preload.data();
            }
            if(bam_iter->ref_id != ref_id || bam_iter->pos >= region_end) {
                if(first_in_region == -1) {
                    return std::vector<uint8_t>();
                }
                return std::vector<uint8_t>(
                        bam_buffer_preload.cbegin() + first_in_region,
                        bam_buffer_preload.cbegin() + (BYTEREF(bam_iter) - bam_buffer_preload.data())
                        );
            }
            bam_iter = BAMREF(bam_next);
            bam_next = BYTEREF(bam_iter) + bam_iter->block_size + 4;
        }
    }

    // last read is still within [region_start, region_end)
    // append data starting from end_intv
    auto ioffset = index.ref[ref_id].ioffset[end_intv];
    auto bgzf_it = bgzf_slicer_iterator_t<SLICER_T>(data, index_coffset(ioffset));
    auto bgzf_end = bgzf_it.end();

    std::vector<uint8_t> bam_buffer;

    bool found_outbound = false;
    size_t next_offset = index_uoffset(ioffset);

    while(bgzf_it != bgzf_end && !found_outbound) {

        bam_iter = reinterpret_cast<const bam_rec_t *>(bam_buffer.data() + next_offset);

        // make sure buffer contains block_size of current read
        size_t min_size = next_offset + 4;
        extend_buffer(bgzf_it, bgzf_end, bam_buffer, min_size);
        bam_iter = reinterpret_cast<const bam_rec_t *>(bam_buffer.data() + next_offset);

        // make sure buffer contains current read
        min_size += bam_iter->block_size + 4;
        extend_buffer(bgzf_it, bgzf_end, bam_buffer, min_size);
        bam_iter = reinterpret_cast<const bam_rec_t *>(bam_buffer.data() + next_offset);

        // try to find out-of-bound read in bam_buffer
        //decltype(bam_iter) bam_prev;
        auto buffer_end = bam_buffer.data() + bam_buffer.size();

        // loop while current read is within buffer
        while(true) {
            if(BYTEREF(bam_iter) + 4 >= buffer_end) break; // read block size is out of bound
            if(BYTEREF(bam_iter) + bam_iter->block_size + 4 >= buffer_end) break; // read end is out of bound

            // read is within buffer
            if(bam_iter->ref_id != ref_id || bam_iter->pos >= region_end) {
                found_outbound = true;
                break;
            }
            if( READ_IN_REGION(bam_iter, region_start, region_end) && first_in_region == -1 ) {
                first_in_region = BYTEREF(bam_iter) - bam_buffer.data() - index_uoffset(ioffset);
            }

            //if(BYTEREF(BAM_NEXT(bam_iter)) >= buffer_end) break; // next read out of bound{
            //    bam_prev = bam_iter;
            //    break;
            //}
            //bam_prev = bam_iter;

            bam_iter = BAM_NEXT(bam_iter);
        }

        // did not find out-of-bound read, readjust prev_offset
        //prev_offset = BYTEREF(bam_prev) - bam_buffer.data();
        next_offset = BYTEREF(bam_iter) - bam_buffer.data(); // calculate next read offset
    }
    
    // concatenate preload buffer with bam_buffer(offset with first uoffset
    bam_buffer_preload.insert(bam_buffer_preload.end(),
            bam_buffer.cbegin() + index_uoffset(ioffset),
            bam_buffer.cbegin() + (BYTEREF(bam_iter) - bam_buffer.data()));

    if(first_in_region == -1)
        return std::vector<uint8_t>();

    return std::vector<uint8_t>(
            bam_buffer_preload.cbegin() + first_in_region,
            bam_buffer_preload.cend());

}


inline uint64_t calc_ioffset(uint64_t coffset, uint64_t uoffset) {
    return (coffset << 16) | uoffset;
}



template<typename SLICER_T>
std::vector<region> bam_to_regions(SLICER_T slicer, size_t start_offset, size_t end_offset) {

    const size_t CHUNK_SZ = 1024*1024;

    auto num_chunks = (end_offset - start_offset) / CHUNK_SZ;
    
    std::vector<size_t> region_starts(num_chunks);

#pragma omp parallel for
    for (size_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {

        auto start_idx = start_offset + (chunk_idx * CHUNK_SZ);


        auto block_idx = bgzf_find_next_block(slicer, start_idx);

        bgzf_slicer_iterator_t<SLICER_T> bgzf_it(slicer, block_idx);
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

#endif
