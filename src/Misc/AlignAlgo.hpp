#ifndef ALIGNALGO_HPP
#define ALIGNALGO_HPP

extern "C" {
    //version lineal
    void* SW_ASM_LIN(Alignment* alignment);
    void* NW_ASM_LIN(Alignment* alignment);
    //version paralelizada SSE
    void* NW_ASM_SSE(Alignment* alignment, bool debug);
    void* SW_ASM_SSE(Alignment* alignment, bool debug);
    //version paralelizada AVX
    void* NW_ASM_AVX(Alignment* alignment, bool debug);
    //void* SW_ASM_AVX(Alignment* alignment, bool debug);
}
#endif  