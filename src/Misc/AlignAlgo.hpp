#ifndef ALIGNALGO_HPP
#define ALIGNALGO_HPP

extern "C" {
    void* SW_ASM_LIN(Alignment* alignment);
    void* NW_ASM_LIN(Alignment* alignment);
    void* NW_ASM_SSE(Alignment* alignment, bool debug);
}
#endif  