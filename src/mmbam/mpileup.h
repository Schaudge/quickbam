#ifndef MMBAM_MPILEUP_H
#define MMBAM_MPILEUP_H

#include <vector>
#include <list>
#include <cstring>
#include <functional>
#include "bam.h"
#include "index.h"
#include <tbb/enumerable_thread_specific.h>

using pileup_filter_func_t = bool (*)(bam_rec_t* read);

//! Read-wise information at a piled up location
struct pileup_info_t {
    size_t buffer_id; /**< which buffer the read was from */
    size_t offset;    /**< where the read is located in the corresponding buffer */
    uint32_t qpos;    /**< the query position of the read that overlaps the piled up position */
    bool is_deletion; /**< whether the read contains a deletion at the piled up position */
};

#define PILEUP_MAX_DEPTH 80000

//! Information about a piled up location, produced by the mpileup engine and offered to the visitor function as a parameter
struct mpileup_t {
    int32_t ref_id;  /**< the numeric reference contig id of the piled up location */
    int32_t pos;     /**< the 0-based genomic coordinate of the piled up location */
    std::vector<const std::vector<uint8_t>*> reads_buffer; /**< The buffers (1 per file) of all reads */
    pileup_info_t *info; /**< additional, read-wise information at the piled up location, length= # of buffers * MAX_PILEUP_DEPTH */
    size_t *depth;       /**< number of reads overlapping this position, length = # of buffers */
    size_t n_buffers;

    //! returns the correct info object
    //! \param b the index of the BAM file (or input buffer)
    //! \param d the index of the read within the specified buffer
    //! \return The info object corresponding to the specified buffer and read
    pileup_info_t& get_info(size_t b, size_t d) {
        return *(info + b * PILEUP_MAX_DEPTH + d);
    }

    //! returns the correct info object (const)
    //! \param b the index of the BAM file (or input buffer)
    //! \param d the index of the read within the specified buffer
    //! \return The info object corresponding to the specified buffer and read
    const pileup_info_t& get_info(size_t b, size_t d) const {
        return *(info + b * PILEUP_MAX_DEPTH + d);
    }

    void allocate(size_t buffers) {
        n_buffers = buffers;
        reads_buffer.resize(n_buffers);
        info = new pileup_info_t[n_buffers * PILEUP_MAX_DEPTH];
        //info = alloc.allocate(n_buffers * PILEUP_MAX_DEPTH);
        depth = new size_t[n_buffers];
        //std::cerr<<"pileup allocated "<<n_buffers<<" n_buffers"<<std::endl;
        //info = tbb::scalable_allocator<pileup_info_t>(n_buffers * PILEUP_MAX_DEPTH);
        //depth = tbb::scalable_allocator<size_t>(n_buffers);
    }

    void set_buffers(const std::vector<std::vector<uint8_t>>& buffers) {
        //std::cerr<<"buffers set"<<std::endl;
        for(size_t i=0; i<n_buffers; i++) reads_buffer[i] = &buffers[i];
    }

    void set_location(int32_t ref_id, int32_t pos, std::vector<std::vector<uint8_t>>& buffers)
    {
        memset(depth, 0, sizeof(size_t) * buffers.size());
        this->ref_id = ref_id;
        this->pos = pos;
    }

    void release() {
        //std::cerr<<"pileup released"<<std::endl;
        if(info != nullptr) delete [] info;
        //alloc.deallocate(info, reads_buffer.size() * PILEUP_MAX_DEPTH);
        if(depth != nullptr) delete [] depth;
        //tbb::scalable_allocator<pileup_info_t>(info, reads_buffer.size() * PILEUP_MAX_DEPTH);
        //tbb::scalable_allocator<size_t>(depth, reads_buffer.size());
    }

};

struct thread_pileup_storage {
    mpileup_t pileups[400];

    thread_pileup_storage(size_t n_buffers) {
        //std::cerr<<"ets allocating pileups"<<std::endl;
        for(size_t i=0; i<400; i++) {
            pileups[i].allocate(n_buffers);
        }
    }

    ~thread_pileup_storage() {
        //std::cerr<<"ets releasing pileups"<<std::endl;
        //for(size_t i=0; i<400; i++) {
        //    pileups[i].release();
        //}
    }

    void set_buffers(const std::vector<std::vector<uint8_t>>& buffers) {
        //std::cerr<<"ets setting buffers"<<std::endl;
        for(size_t i=0; i<400; i++) {
            pileups[i].set_buffers(buffers);
        }
    }
};

tbb::enumerable_thread_specific<thread_pileup_storage> mpileup_ets() {
    return tbb::enumerable_thread_specific<thread_pileup_storage>(2);
}

//! Type alias of multiple input mfile
using mfiles_t = std::vector<std::reference_wrapper<const mfile_t::ptr_t>>;

//! Type alias of multiple index
using indices_t = std::vector<std::reference_wrapper<const index_t>>;

constexpr size_t max_buff = 256 * 1024 * 1024;

const bam_rec_t* multi_buffer_next_read(
        std::vector<std::vector<uint8_t>>& m_buffer,
        std::vector<size_t>& m_buffer_offsets,
        size_t& buffer_id) {
    size_t which_buffer = 0;
    int32_t smallest_pos = 0;
    const bam_rec_t* first_read = nullptr;
    for(size_t i_buffer = 0; i_buffer<m_buffer.size(); i_buffer++) {
        auto offset = m_buffer_offsets[i_buffer];
        if(offset < m_buffer[i_buffer].size()) {
            // i_buffer has read at offset
            const bam_rec_t *read = BAMREF(m_buffer[i_buffer].data() + offset);
            if(first_read == nullptr || read->pos < smallest_pos) {
                which_buffer = i_buffer;
                smallest_pos = read->pos;
                first_read = read;
            }
        }
    }

    if(first_read != nullptr) {
        // read found, increment buffer offset
        m_buffer_offsets[which_buffer] += (first_read->block_size + 4);
    }

    buffer_id = which_buffer;


    return first_read;
}

mpileup_t* allocate_mpileups(const std::vector<std::vector<uint8_t>>& m_buffers) {
    mpileup_t *pileups = new mpileup_t[400];
    for(size_t i=0; i<400; i++) {
        pileups[i].allocate(m_buffers.size());
        pileups[i].set_buffers(m_buffers);

    }
    return pileups;
}

void release_mpileups(mpileup_t* pileups) {
    for(size_t i=0; i<400; i++) pileups[i].release();
    delete [] pileups;
}

//! Multiple input pileup engine
//
//! \param mfiles A list of input memory mapped BAM files
//! \param indices A list of indices of the input BAM files, order-matched
//! \param ref_id The numerical reference contig id of the pileup region
//! \param pos_start The genomic coordinate of the start of the pileup region
//! \param pos_end The genomic coordinate of the end of the pileup region
//! \param predicate_func A functor, invoked with one parameter of a bam_rec_t struct const reference, and returns true if such read should be considered by the pileup engine
//! \param visitor_func A functor, invoked at each piled up location with one parameter of a mpileup_t struct const reference, and returns true if pileup engine should continue


template <typename SLICER_T>
void mpileup(std::vector<std::reference_wrapper<const mfile_t::ptr_t>> mfiles,
        std::vector<SLICER_T> slicers,
        std::vector<std::reference_wrapper<const index_t>> indices,
        tbb::enumerable_thread_specific<thread_pileup_storage>& ets,
        uint32_t ref_id, int32_t pos_start, int32_t pos_end,
        const std::function<bool(const bam_rec_t&)>& predicate_func,
        const std::function<bool(const mpileup_t&)>& visitor_func) {

    const auto n_files = mfiles.size();

    if(n_files != indices.size()) throw std::runtime_error("number of files and indices mismatch");

    // load reads in region from all files
    std::vector<std::vector<uint8_t>> m_buffers(n_files);
    std::vector<size_t> m_buffer_offsets(n_files, 0);
    for(size_t i=0; i<n_files; i++) {
        //mfile_slicer_t data(mfiles[i]);
        m_buffers[i] = bam_load_region(slicers[i], mfiles[i], indices[i], ref_id, pos_start, pos_end);
    }

    // initialize pileup data structure
    bool exists;
    auto this_ets = ets.local(exists);
    //std::cerr<<"have existing ets object? "<<exists<<std::endl;
    //std::cerr<<"acquired thread local storage object"<<std::endl;
    //auto pileups = allocate_mpileups(m_buffers);
    this_ets.set_buffers(m_buffers);
    //std::cerr<<"set thread local storage buffers"<<std::endl;
    auto * pileups = &(this_ets.pileups[0]);
    int p_head = 0;
    int p_tail = 0;
    int p_size = 400;

    // pileup until reads are exhausted
    size_t buffer_id;
    const bam_rec_t *bam_it;
    while((bam_it = multi_buffer_next_read(m_buffers, m_buffer_offsets, buffer_id)) != nullptr) {
        if(!READ_IN_REGION(bam_it, pos_start, pos_end)) { continue; }
        if(!predicate_func(*bam_it)) { continue; }

        while(p_head != p_tail) {
            if(pileups[p_head].ref_id < bam_it->ref_id || pileups[p_head].pos < bam_it->pos) {
                if(!visitor_func(pileups[p_head])) {
                    //release_mpileups(pileups);
                    return;
                }
                p_head++;
                if(p_head == p_size) p_head = 0;
            }
            else break;
        }

        int p_iter = p_head;

        // walk through the read, and update pileups
        const uint32_t *cigar_ops = (const uint32_t *)(
                (const uint8_t *)(bam_it)
                + sizeof(bam_rec_t) 
                + sizeof(char) * bam_it->l_read_name );

        const uint8_t *bam_begin = m_buffers[buffer_id].data();
        size_t read_offset = BYTEREF(bam_it) - bam_begin;
        uint32_t qpos = 0;
        int32_t rpos = bam_it->pos;

        for(int op=0; op < bam_it->n_cigar_op; op++) {
            auto cigar_len = cigar_ops[op] >> 4;
            auto cigar_op  = cigar_ops[op] & 0xF;
            switch(cigar_op) {
                case 1: // I: insertion only consumes query 
                case 4: // S: soft padded only consumes query
                    qpos += cigar_len;
                    break;
                case 2: // D: deletion only consumes reference
                case 3: // N: reference skip only consumes reference
                    while(cigar_len > 0) {
        
                        if(p_iter != p_tail) {
                            mpileup_t *p = &pileups[p_iter];
                            if(p->depth[buffer_id] < PILEUP_MAX_DEPTH) p->get_info(buffer_id, p->depth[buffer_id]++) = {
                                .buffer_id = buffer_id,
                                .offset = read_offset, 
                                .qpos = qpos,
                                .is_deletion = true};
                            p_iter++;
                            if(p_iter == p_size) p_iter = 0;
                        }
                        else {
                            int next_tail = p_tail + 1;
                            if(next_tail == p_size) next_tail = 0;
                            if(next_tail != p_head) {
                                mpileup_t *p = &pileups[p_tail];
                                p->set_location(bam_it->ref_id, rpos, m_buffers);
                                p->get_info(buffer_id, p->depth[buffer_id]++) = {
                                    .buffer_id = buffer_id,
                                    .offset = read_offset,
                                    .qpos = qpos,
                                    .is_deletion = true};
                                p_tail = next_tail;
                                p_iter = p_tail;
                            } else {
                                std::cerr<<"circular buffer full, result truncated"<<std::endl;
                            }
                        }

                        rpos++;
                        cigar_len--;
                    }
                    break;
                case 0: // M: alignment match consumes both
                case 7: // =: sequence match consumes both
                case 8: // X: sequence mismatch consumes both

                    while(cigar_len > 0) {
        
                        if(p_iter != p_tail) {
                            mpileup_t *p = &pileups[p_iter];
                            if(p->depth[buffer_id] < PILEUP_MAX_DEPTH) {
                                /*
                                p->get_info(buffer_id, p->depth[buffer_id]++) = {
                                .buffer_id = buffer_id,
                                .offset = read_offset, 
                                .qpos = qpos,
                                .is_deletion = false};
                                */
                                auto& info = p->get_info(buffer_id, p->depth[buffer_id]++);
                                info.buffer_id = buffer_id;
                                info.offset = read_offset;
                                info.qpos = qpos;
                                info.is_deletion = false;
                            }
                            p_iter++;
                            if(p_iter == p_size) p_iter = 0;
                        }
                        else {
                            int next_tail = p_tail + 1;
                            if(next_tail == p_size) next_tail = 0;
                            if(next_tail != p_head) {
                                mpileup_t *p = &pileups[p_tail];
                                p->set_location(bam_it->ref_id, rpos, m_buffers);
                                p->get_info(buffer_id, p->depth[buffer_id]++) = {
                                    .buffer_id = buffer_id,
                                    .offset = read_offset,
                                    .qpos = qpos,
                                    .is_deletion = false};
                                p_tail = next_tail;
                                p_iter = p_tail;
                            }
                            else {
                                std::cerr<<"circular buffer full, result truncated"<<std::endl;
                            }
                        }

                        qpos++;
                        rpos++;
                        cigar_len--;
                    }
                    break;
                default: // other operations consume neither
                    break;
            }
        }
    }

    while(p_head != p_tail) {
        // emit the rest
        if(p_head == p_size) p_head = 0;
        if(p_head == p_tail) break;

        if(pileups[p_head].pos >= pos_end) {
            //release_mpileups(pileups);
            return;
        }

        if(!visitor_func(pileups[p_head++])) {
            //release_mpileups(pileups);
            return;
        }
    }

    //release_mpileups(pileups);
}
#endif
