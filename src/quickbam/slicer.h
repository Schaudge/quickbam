#ifndef QUICKBAM_SLICER_H
#define QUICKBAM_SLICER_H

#include <memory>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <sys/stat.h>

#include <quickbam/mfile.h>

/*! Slicer type that abstracts over an mfile */
struct mfile_slicer_t {
    using ptr_t = std::shared_ptr<const uint8_t[]>;
    const mfile_t::ptr_t& mfile;

    mfile_slicer_t(const mfile_t::ptr_t& mfile) : mfile(mfile)
    {
    }

    ptr_t slice(uint64_t start, uint64_t) {
        auto buf = begin<const uint8_t>(mfile) + start;
        ptr_t ptr(buf, mfile_deleter());
        // TODO: Isn't this move implicit?
        return std::move(ptr);
    }

    size_t size() {
        return mfile->size;
    }

    bool operator==(const mfile_slicer_t& rhs) { return mfile == rhs.mfile; }
};



//! A slicer that abstracts over native filesystem operations (pread)
/*! Since many BAM files are too large to fit in memory, quickbam uses an
 * abstraction called "slicers". These are simply types that provide the
 * 'slice' and 'size' methods. The purpose of this abstraction is to allow an
 * entire BAM file to be operated on, while provided the owner the ability to
 * control which data is in memory at any given time.
 *
 * The file_slicer_t implements the slicer interface for generic files
 * accessed directly through the filesystem, including systems like NFS and
 * Lustre.
 */
struct file_slicer_t {

    //! type alias to the smart pointer type
    using ptr_t = std::shared_ptr<const uint8_t[]>;

    size_t file_size;
    int fd;

    //! Creates a file slicer from a string path
    file_slicer_t(std::string file_path)
    {
        struct stat stat_;
        auto rc = stat(file_path.c_str(), &stat_);
        file_size = stat_.st_size;

        fd = open(file_path.c_str(), O_RDONLY);

        if (fd == -1) {
            // TODO: We might want to consider changing the interface of
            // file_slicer_t to take fd directly, and creating a wrapper
            // function that can fail without throwing an exception.
            throw std::invalid_argument("Failed to open file: '" + file_path +
                "'. Are you sure the path is correct?");
        }
    }

    //! Returns a slice containing the requested byte range of the underlying file
    /*! \param start the 0-indexed byte offset the returned slice should start at
     *  \param end the 0-indexed byte offset the returned slice should end at (final byte is not returned)
     *  \return a ptr_t slice of bytes
     */
    ptr_t slice(uint64_t start, uint64_t end) {
        size_t range_size = end - start;
        ptr_t buf(new uint8_t[range_size]);
        uint8_t *buf_ptr = (uint8_t*)buf.get();

        while(range_size > 0) {
            auto n_read = pread(fd, buf_ptr, range_size, start);
            if (n_read == 0) {
                break;
            }
            range_size -= n_read;
            start += n_read;
            buf_ptr += n_read;
        }

        return buf;
    }

    //! Returns the size of the entire underlying file
    /*! \return the size of the underlying file in bytes
     */
    size_t size() {
        return file_size;
    }

    bool operator==(const file_slicer_t& rhs) { return fd == rhs.fd && file_size == rhs.file_size; }
};

#endif // SLICER_H

