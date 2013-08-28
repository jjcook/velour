//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// bufferedio.h
//

#ifndef VELOUR_BUFFEREDIO_H
#define VELOUR_BUFFEREDIO_H

class AlignedWriteBuffer
{
    static const size_t buffer_size = 16384; // XXX: don't change this, as its hard coded in other code!
    static const size_t n_buffer = 128; // 2MB buffer

    char buffer_[n_buffer][buffer_size];
    size_t next_buffer_;
    size_t buffer_offset_;

    int filedes_;

    public:
        AlignedWriteBuffer(char *filename);
        ~AlignedWriteBuffer();

        void * RequestBuffer(size_t amount);
};

AlignedWriteBuffer::AlignedWriteBuffer(char *filename) : next_buffer_(0), buffer_offset_(0), filedes_(0)
{
    filedes_ = open(filename, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);

    if (filedes_ == -1)
    {
        fprintf(stderr, "Failed to open output file: %s\n", filename);
        perror("reason: ");
        exit(EXIT_FAILURE);
    }
}

AlignedWriteBuffer::~AlignedWriteBuffer()
{
    if (buffer_offset_ < buffer_size) {
        // XXX: record wasted space
        buffer_[next_buffer_][buffer_offset_] = -1; // end-of-buffer delimiter if space left
    }
    ++ next_buffer_;

    size_t write_bytes = next_buffer_ * buffer_size;
    ssize_t retval = write(filedes_, static_cast<void*>(buffer_), write_bytes);
    if (retval == -1) {
        fprintf(stderr, "Failed to write to output file.\n");
        perror("reason: ");
        exit(EXIT_FAILURE);
    }

    close(filedes_);
}
        
void * AlignedWriteBuffer::RequestBuffer(size_t amount)
{
    if (amount > buffer_size) {
        fprintf(stderr, "INTERNAL ERROR: Buffer not large enough for requested quanta: %zu in %zu bytes.\n", amount, buffer_size);
        exit(EXIT_FAILURE);
    }

    if (buffer_offset_ + amount > buffer_size) { // amount to write overruns current aligned buffer
        if (buffer_offset_ < buffer_size) {
            // XXX: record wasted space
            buffer_[next_buffer_][buffer_offset_] = -1; // end-of-buffer delimiter if space left
        }

        ++ next_buffer_;
        buffer_offset_ = 0;
    }

    if (next_buffer_ == n_buffer) {
        size_t write_bytes = n_buffer * buffer_size;
        ssize_t retval = write(filedes_, static_cast<void*>(buffer_), write_bytes);
        if (retval == -1) {
            fprintf(stderr, "Failed to write to output file.\n");
            perror("reason: ");
            exit(EXIT_FAILURE);
        }

        // then, reset
        next_buffer_ = 0;
    }

    void * retval = &buffer_[next_buffer_][buffer_offset_];

    buffer_offset_ += amount;

    return retval;
}

#endif // VELOUR_BUFFEREDIO_H
