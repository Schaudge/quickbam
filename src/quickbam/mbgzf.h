#ifndef QUICKBAM_MBGZF_H
#define QUICKBAM_MBGZF_H

/****************************************************************
 * BGZF APIs implemented on memory buffers
*****************************************************************/
#include <string.h>
#include <vector>
#include <iterator>
#include <functional>
#include "mfile.h"
#include "nfo_iterator.h"
#include <cassert>

/****************************************************************
 * API usage:
 *
 * To create a bgzf proxy over an mfile
 *      auto mfile = mfile_open("filename");
 *      bgzf_mfile_proxy_t bgzf_proxy(mfile);
 *
 * To create bgzf block iterators from a bgzf proxy
 *      auto it_begin = bgzf_proxy.begin();
 *      auto it_end   = bgzf_proxy.end();
 *
 * To check if a bgzf block is the EOF marker block
 *      auto is_eof = is_bgzf_eof_block(*it_begin)
 *
 * To iterate over an entire file using range for
 *      for(auto& block : bgzf_proxy) {
 *          // e.g. inflate the block
 *          auto inflated = bgzf_inflate(block)
 *      }
*****************************************************************/

constexpr auto bgzf_xlen_offset = 10;
constexpr auto bgzf_block_size_offset = 16;
constexpr auto bgzf_block_size_size = sizeof(uint16_t);
constexpr auto bgzf_cdata_offset = 18;

const uint32_t BGZF_MAX_BLOCK_SIZE = 65536;

constexpr uint8_t eof_bytes[] = {
    0x1f, 0x8b, 0x08, 0x04,
    0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00,
    0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00,
    0,0,0,0,0,0,0,0
};

//! This struct mirrors the byte layout of the BGZF blocks inside BAM files, per SAM file specification section 4.1
struct bgzf_block_t {
    uint8_t id1; uint8_t id2; uint8_t cm; uint8_t flg;
    uint32_t mtime; uint8_t xfl; uint8_t os; uint16_t xlen;
    uint8_t si1; uint8_t si2; uint16_t slen; uint16_t bsize;
    uint8_t cdata[];

    // uint32_t crc32;
    // uint32_t isize;

    bgzf_block_t(const bgzf_block_t&) = delete;
};

//! checks if a bgzf block is the special EOF block
//! \param block The bgzf block to check
//! \return True if the bgzf block is the special EOF block, false otherwise
inline bool is_bgzf_eof_block(const bgzf_block_t& block) {
    return strncmp((const char*)eof_bytes, (const char*)&block, 28) == 0;
}

//! BGZF iterator, which is defined as a specialization of nfo_iterator_t
using bgzf_iterator_t = nfo_iterator<bgzf_block_t, uint16_t, bgzf_block_size_offset, 1>;

//! Iterator which takes its input from a slicer and iterates over bgzf_block_t elements.
template<typename SLICER_T>
struct bgzf_slicer_iterator_t : std::iterator<std::forward_iterator_tag, bgzf_block_t> {
    SLICER_T data;
    typename SLICER_T::ptr_t current_block_slice;
    size_t current_offset;
    size_t current_block_size;

    void load_current_block() {
        if ((current_offset + bgzf_block_size_offset + bgzf_block_size_size) >= data.size()) {
            return;
        }

        auto block_size_slice_start = current_offset + bgzf_block_size_offset;
        auto block_size_slice_end = block_size_slice_start + bgzf_block_size_size;
        auto block_size_slice = data.slice(block_size_slice_start, block_size_slice_end);
        auto block_size = *reinterpret_cast<const uint16_t*>(block_size_slice.get()) + 1;

        current_block_size = block_size;

        auto block_slice_start = current_offset;
        auto block_slice_end = current_offset + block_size;
        current_block_slice = data.slice(block_slice_start, block_slice_end);

        // Special case for iterator end()
        if (current_offset != data.size()) {
            assert(current_block_slice[0] == 31);
            assert(current_block_slice[1] == 139);
        }
    }

    void load_next_block() {
        current_offset += current_block_size;
        current_block_size = 0;
        load_current_block();
    }


    //! Creates iterator from beginning of slicer data
    bgzf_slicer_iterator_t(SLICER_T slicer) : current_offset(0), current_block_size(0), data(slicer) {
        load_current_block();
    }

    //! Create iterator from specific offset of slicer data
    /*! \param slicer Data slicer
     *  \param starting_offset 0-indexed byte offset to start iterating from
     */
    bgzf_slicer_iterator_t(SLICER_T slicer, size_t starting_offset) :
            current_offset(starting_offset), current_block_size(0), data(slicer) {
        load_current_block();
    }

    bgzf_slicer_iterator_t operator++() {
        load_next_block();
        return *this;
    }
    bgzf_slicer_iterator_t operator++(int) {
        load_next_block();
        return *this;
    }
    const bgzf_block_t& operator*() { return *reinterpret_cast<const bgzf_block_t*>(current_block_slice.get()); }

    bool operator==(const bgzf_slicer_iterator_t& rhs) {
        return data == rhs.data && current_offset == rhs.current_offset;
    }
    bool operator!=(const bgzf_slicer_iterator_t& rhs) { return !(*this == rhs); }

    //! Returns iterator at beginning of slicer data
    bgzf_slicer_iterator_t begin() {
        return bgzf_slicer_iterator_t(data, 0);
    }

    //! Returns iterator at end of slicer data
    bgzf_slicer_iterator_t end() {
        return bgzf_slicer_iterator_t(data, data.size());
    }

};

//! Creates a bgzf iterator with the given mfile, and with the record starting at the specified offset
//! \param mfile The mfile from which the iterator should be created
//! \param offset The offset at which the iterator should point at. It is
//important that the offset is an actual record start
//! \return A BGZF iterator
inline auto bgzf_iterator_at(const mfile_t::ptr_t& mfile, off_t offset) {
    return bgzf_iterator_t(reinterpret_cast<const bgzf_block_t *>(::begin(mfile) + offset));
}

//! calculates the inflated data size of the given BGZF block
//! \param block The BGZF block with which to calculate the inflate size
//! \return Number of bytes the inflated data will be.
inline auto bgzf_isize(const bgzf_block_t& block) {
    return *(uint32_t *)((uint8_t*)(&block) + (block.bsize +1) - 4);
}

//! This struct provides methods to create various BGZF iterators from an underlying mfile
struct bgzf_mfile_proxy_t {
    const mfile_t::ptr_t& mfile;

    //! Initializes a bgzf proxy from the given mfile
    //! \param mfile On which file the proxy should be created
    bgzf_mfile_proxy_t(const mfile_t::ptr_t& mfile) : mfile(mfile) {}

    //! Returns a BGZF iterator at the beginning of the file
    bgzf_iterator_t begin() { return bgzf_iterator_t(::begin<bgzf_block_t>(mfile)); }
    //! Returns a BGZF iterator at an arbitary offset of the file
    bgzf_iterator_t at_offset(size_t offset) { return bgzf_iterator_t(reinterpret_cast<const bgzf_block_t *>(::begin(mfile) + offset)); }
    //! Returns a BGZF iterator past the last byte of the file
    bgzf_iterator_t end() { return bgzf_iterator_t(::end<bgzf_block_t>(mfile)); }
};

//! Inflate the given BGZF block
//! \param block The BGZF block to inflate
//! \return A vector consists of the inflated bytes
std::vector<uint8_t> bgzf_inflate(const bgzf_block_t& block);

//! Inflate a continuous range of several BGZF blocks
//! \param src The beginning of the BGZF blocks
//! \param src_len The length of all the BGZF blocks combined
//! \return A vector consists of the inflated bytes from all specified BGZF blocks
std::vector<uint8_t> bgzf_inflate_range(const uint8_t *src, const size_t src_len);

//! Inflate a continuous range of several BGZF blocks in parallel
//! \param src The beginning of the BGZF blocks
//! \param src_len The length of all the BGZF blocks combined
//! \param parallel_inflate A functor object to implement the actual parallelization mechanism.
//! \return A vector consists of the inflated bytes from all specified BGZF blocks
std::vector<uint8_t> bgzf_inflate_range_p(const uint8_t *src, const size_t src_len,
        const std::function<std::vector<uint8_t>(
            const uint8_t*,
            const size_t,
            const size_t,
            const std::vector<size_t>&,
            const std::vector<size_t>&,
            size_t (*)(const uint8_t*, size_t, uint8_t*, size_t))>& parallel_inflate);


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

#endif
