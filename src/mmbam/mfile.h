#ifndef MMAP_FILE_H
#define MMAP_FILE_H

#include <memory>
#include <stdexcept>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <iostream>



/*! Memory mapped file managed by C++ smart-pointer */
struct mfile_t {
    int fd;         /**< file descriptor of the underlying file */
    off_t size;     /**< size of the underlying file */
    uint8_t *mmptr; /**< byte pointer to the memory mapped virtual address */

    //! type alias to the smart pointer type
    using ptr_t = std::unique_ptr<mfile_t>;

    //! constructor to create a memory mapped file from a string
    /*! This constructor will open the file given in the parameter, finds its size,
     *  and memory map it into the program's virtual memory space.
     *  \param filename the name of the file to open
     *  \return a mfile_t struct initialized from opeing the given file via mmap
     *  \throw std::runtime_error if somehow the call to mmap fails
     */
    mfile_t(const char* filename) :
        fd(open(filename, O_RDONLY)), mmptr(nullptr)
    {
        if(fd == -1) throw std::runtime_error("cannot open file");

        // get file size
        struct stat s;
        fstat(fd, &s);
        size = s.st_size;

        mmptr = (uint8_t*)mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0);
        if(mmptr == MAP_FAILED) {
            close(fd);
            throw std::runtime_error("cannot memory map");
        }
    }
    ~mfile_t() { munmap(mmptr, size); close(fd); }

    // delete copy constructors
    mfile_t(const mfile_t&) = delete;
    mfile_t& operator=(const mfile_t&) = delete;
};

//! a convinent function to create a mfile_t from a stringly type
//! \tparam T The type of the underlying records of file
//! \param filename an object of any stringly type that can be implicitly cast to ``const char *``
//! \return a smart unique pointer to mfile_t struct initialized from opening the given file.
template <typename T>
auto mfile_open(T&& filename) {
    return std::make_unique<mfile_t>(std::forward<T>(filename));
}

//! The main access function, which returns a pointer of an arbitary type to act as an in-memory buffer
//! \tparam T The type of the underlying records of file
//! \param mfile The mfile from which a pointer is to be created
//! \return A pointer of an arbitary type, pointed at the beginning of the file
template <typename T = char>
const T* begin(const mfile_t::ptr_t& mfile) noexcept {
    if(mfile->mmptr == nullptr) return nullptr;
    return reinterpret_cast<T *>(mfile->mmptr);
}

//! This function returns a pointer of an arbitary type, which points at the end position of the mfile
//! This is useful together with begin() in a loop to iterate over records in the mfile
//! \tparam T The type of the underlying records of file
//! \param mfile The mfile from which the end-pointer is to be created
//! \return A pointer of an arbitary type, pointed at the end of the file (exactly 1 byte past the last byte)
template <typename T = char>
const T* end(const mfile_t::ptr_t& mfile) noexcept {
    return reinterpret_cast<T*>(mfile->mmptr + mfile->size);
}

// no-op deleter to avoid unique_ptr trying to free the mmap pointer
struct mfile_deleter{  
    void operator()(const uint8_t[]) {
    }
};


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

struct file_slicer_t {

    using ptr_t = std::shared_ptr<const uint8_t[]>;

    size_t file_size;
    int fd;

    file_slicer_t(std::string file_path)
    {
        struct stat stat_;
        auto rc = stat(file_path.c_str(), &stat_);
        file_size = stat_.st_size;

        fd = open(file_path.c_str(), O_RDONLY);
    }

    ptr_t slice(uint64_t start, uint64_t end) {
        const size_t range_size = end - start;
        ptr_t buf(new uint8_t[range_size]);
        auto n_read = pread(fd, (void*)buf.get(), range_size, start);
        // TODO: Isn't this move implicit?
        return std::move(buf);
    }

    size_t size() {
        return file_size;
    }

    bool operator==(const file_slicer_t& rhs) { return fd == rhs.fd && file_size == rhs.file_size; }
};



#endif
