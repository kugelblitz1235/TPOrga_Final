#ifndef ALIGNALGO_HPP
#define ALIGNALGO_HPP

extern "C" {
    void* SWLIN(
        char* seq1,
        char* seq2,
        short match,
        short missmatch,
        short gap,
        short len1,
        short len2,
        Alignment* alignment);
}
#endif  