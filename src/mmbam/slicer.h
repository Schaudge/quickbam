#ifndef SLICER_H
#define SLICER_H

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

/*! Slicer type that abstracts over a file */
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

#endif // SLICER_H

