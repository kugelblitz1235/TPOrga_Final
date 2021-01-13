#include <unistd.h>
#include <iostream>
#include <string>
#include <emmintrin.h>

#include "Misc/littletest.hpp"
#include "Misc/Types.hpp"
#include "Needleman_Wunsch/NW_C.hpp"
#include "Smith_Waterman/SW_C.hpp"
#include "Misc/JSON.hpp"
#include "Misc/AlignAlgo.hpp"
#include "Misc/Utility.hpp"


#define DBG(x) cout << #x << " = " << (x) <<"\n"
using namespace std;

char* random_seq(const int len) {

    char *s =(char *)malloc(len * sizeof(char) + 1);
  
    static const char bases[] = "ACTG";

    for (int i = 0; i < len; ++i) 
        s[i] = bases[rand() % (sizeof(bases) - 1)];
    
    s[len] = '\0';

    return s;
}

//===========================Needleman wunsch===========================================
LT_BEGIN_SUITE(TestNW)

void set_up() {
    printf("Ejecutando test para NW...\n");
}
void tear_down() {}

LT_END_SUITE(TestNW)

LT_BEGIN_TEST(TestNW, NW_C_x5)
    int vector_len = 32;
    int s_len1 = 32 + rand() % 1000;
    int s_len2 = 32 + rand() % 1000;

    char *s1 = random_seq(s_len1);
    char *s2 = random_seq(s_len2);

    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);
    Alignment* alignment_logic = new_alignment();
    alignment_logic->sequence_1 = new_Sequence_from_string(s1);
    alignment_logic->sequence_2 = new_Sequence_from_string(s2);
    Alignment* alignment_sse = new_alignment();
    alignment_sse->sequence_1 = new_Sequence_from_string(s1);
    alignment_sse->sequence_2 = new_Sequence_from_string(s2);
    Alignment* alignment_avx = new_alignment();
    alignment_avx->sequence_1 = new_Sequence_from_string(s1);
    alignment_avx->sequence_2 = new_Sequence_from_string(s2);
    Alignment* alignment_avx512 = new_alignment();
    alignment_avx512->sequence_1 = new_Sequence_from_string(s1);
    alignment_avx512->sequence_2 = new_Sequence_from_string(s2);

    // Actualizo valor de s_len para considerar el -
    s_len1 = alignment_sse->sequence_1->length;
    s_len2 = alignment_sse->sequence_2->length;

    NW::C::LIN::NW(*alignment_lin, true);
    NW::C::with_logic_SSE(*alignment_logic, vector_len, true);
    NW::C::SSE::NW(*alignment_sse, true);
    NW::C::AVX::NW(*alignment_avx, true);
    NW::C::AVX512::NW(*alignment_avx512, true);

    bool score = (alignment_avx512->result->score == alignment_lin->result->score) &&
                 (alignment_avx->result->score == alignment_lin->result->score) &&
                 (alignment_sse->result->score == alignment_lin->result->score) &&
                 (alignment_logic->result->score == alignment_lin->result->score);

    LT_CHECK(score);

    bool valid_seqs = valid_alignment(*alignment_lin) &&
                      valid_alignment(*alignment_logic) &&
                      valid_alignment(*alignment_sse) &&
                      valid_alignment(*alignment_avx) &&
                      valid_alignment(*alignment_avx512);

    LT_CHECK(valid_seqs);

    bool matrix_ok = true;
    bool position_ok;
    for(int i = 0 ; i < s_len2 && matrix_ok ; i++){
        for(int j = 0 ; j < s_len1 && matrix_ok ; j++){
            position_ok = get_score_SSE(alignment_sse->matrix,s_len1,i,j,8) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8) &&
                          get_score_SSE(alignment_logic->matrix,s_len1,i,j,32) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8) && 
                          get_score_SSE(alignment_avx->matrix,s_len1,i,j,16) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8) && 
                          get_score_SSE(alignment_avx512->matrix,s_len1,i,j,32) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8);
            matrix_ok &= position_ok;
            if (!position_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_LIN(alignment_lin->matrix,s_len1,i,j,16));
                DBG(get_score_SSE(alignment_logic->matrix,s_len1,i,j,32));
                DBG(get_score_SSE(alignment_sse->matrix,s_len1,i,j,8));
                DBG(get_score_SSE(alignment_avx->matrix,s_len1,i,j,16));   
                DBG(get_score_SSE(alignment_avx512->matrix,s_len1,i,j,32));  
            }
        }
    }
    
    LT_CHECK( matrix_ok );

    if (!matrix_ok || !valid_seqs || !score) {
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score Logic: %d\n", alignment_logic->result->score);
        printf("Score SSE: %d\n", alignment_sse->result->score);
        printf("Score AVX: %d\n", alignment_avx->result->score);
        printf("Score AVX512: %d\n", alignment_avx512->result->score);

        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment Logic:\n");
        printf("%s\n", alignment_logic->result->sequence_1->sequence);
        printf("%s\n", alignment_logic->result->sequence_2->sequence);

        printf("Alignment SSE:\n");
        printf("%s\n", alignment_sse->result->sequence_1->sequence);
        printf("%s\n", alignment_sse->result->sequence_2->sequence);

        printf("Alignment AVX:\n");
        printf("%s\n", alignment_avx->result->sequence_1->sequence);
        printf("%s\n", alignment_avx->result->sequence_2->sequence);
        
        printf("Alignment AVX512:\n");
        printf("%s\n", alignment_avx512->result->sequence_1->sequence);
        printf("%s\n", alignment_avx512->result->sequence_2->sequence);
        
        ofstream ofs_logic("NW_C_score_matrix_logic.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_logic->matrix, alignment_logic, 32, ofs_logic);
        ofs_logic.close();
        ofstream ofs_SSE("NW_C_score_matrix_SSE.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_sse->matrix, alignment_sse, 8, ofs_SSE);
        ofs_SSE.close();
        ofstream ofs_AVX("NW_C_score_matrix_AVX.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_avx->matrix, alignment_avx, 16, ofs_AVX);
        ofs_AVX.close();
        ofstream ofs_AVX512("NW_C_score_matrix_AVX512.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_avx512->matrix, alignment_avx512, 32, ofs_AVX512);
        ofs_AVX512.close();
    }

    free(alignment_avx512->matrix);
    free(alignment_avx->matrix);
    free(alignment_sse->matrix);
    free(alignment_logic->matrix);
    free(alignment_lin->matrix);
LT_END_TEST(NW_C_x5)

LT_BEGIN_TEST(TestNW, NW_ASM_SSE_test)
    int vector_len = 8;
    int s_len1 = vector_len + rand() % 100;
    int s_len2 = vector_len + rand() % 100;
    char *s1 = random_seq(s_len1);
    char *s2 = random_seq(s_len2);

    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string(s1);
    alignment->sequence_2 = new_Sequence_from_string(s2);

    // Actualizo valor de s_len para considerar el -
    s_len1 = alignment->sequence_1->length;
    s_len2 = alignment->sequence_2->length;

    NW::C::LIN::NW(*alignment_lin, true);
    
    NW_ASM_SSE(alignment, true);

    bool score = (alignment->result->score == alignment_lin->result->score);
    LT_CHECK(score);

    bool valid_seqs = valid_alignment(*alignment_lin) && valid_alignment(*alignment);

    LT_CHECK(valid_seqs);

    bool matrix_ok = true;
    bool position_ok;
    for(int i = 0 ; i < s_len2 ; i++){
        for(int j = 0 ; j < s_len1 ; j++){
            position_ok = get_score_SSE(alignment->matrix,s_len1,i,j,vector_len) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,vector_len);
            matrix_ok &= position_ok;
            if (!position_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_LIN(alignment_lin->matrix,s_len1,i,j,vector_len));
                DBG(get_score_SSE(alignment->matrix,s_len1,i,j,vector_len));   
            }
        }
    }

 
    LT_CHECK( matrix_ok );

    if (!matrix_ok || !valid_seqs || !score) {
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score SSE: %d\n", alignment->result->score);

        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment SSE:\n");
        printf("%s\n", alignment->result->sequence_1->sequence);
        printf("%s\n", alignment->result->sequence_2->sequence);

        // ofstream ofs_C_SSE("NW_C_score_matrix_SSE.txt", std::ofstream::trunc);
        // printScoreMatrix(alignment_lin->matrix, alignment_lin, 8, ofs_C_SSE);
        // ofs_C_SSE.close();
        ofstream ofs_SSE("NW_ASM_score_matrix_SSE.txt", std::ofstream::trunc);
        printScoreMatrix(alignment->matrix, alignment, vector_len, ofs_SSE);
        ofs_SSE.close();
    }

    free(alignment->matrix);
    free(alignment_lin->matrix);
LT_END_TEST(NW_ASM_SSE_test)

LT_BEGIN_TEST(TestNW, NW_ASM_AVX_test)
    int vector_len = 16;
    int s_len1 = vector_len + rand() % 100;
    int s_len2 = vector_len + rand() % 100;
    
    char *s1 = random_seq(s_len1);
    char *s2 = random_seq(s_len2);

    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string(s1);
    alignment->sequence_2 = new_Sequence_from_string(s2);

    // Actualizo valor de s_len para considerar el -
    s_len1 = alignment->sequence_1->length;
    s_len2 = alignment->sequence_2->length;

    NW::C::LIN::NW(*alignment_lin, true);
    
    NW_ASM_AVX(alignment, true);

    bool score = (alignment->result->score == alignment_lin->result->score);
    LT_CHECK(score);

    bool valid_seqs = valid_alignment(*alignment_lin) && valid_alignment(*alignment);

    LT_CHECK(valid_seqs);

    bool matrix_ok = true;
    bool position_ok;
    for(int i = 0 ; i < s_len2 ; i++){
        for(int j = 0 ; j < s_len1 ; j++){
            position_ok = get_score_SSE(alignment->matrix,s_len1,i,j,vector_len) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,vector_len);
            matrix_ok &= position_ok;
            if (!position_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_LIN(alignment_lin->matrix,s_len1,i,j,8));
                DBG(get_score_SSE(alignment->matrix,s_len1,i,j,vector_len));   
            }
        }
    }
    
    LT_CHECK( matrix_ok );

    if (!matrix_ok || !valid_seqs || !score) {
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score AVX: %d\n", alignment->result->score);

        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment AVX:\n");
        printf("%s\n", alignment->result->sequence_1->sequence);
        printf("%s\n", alignment->result->sequence_2->sequence);

        ofstream ofs_AVX("NW_ASM_score_matrix_AVX.txt", std::ofstream::trunc);
        printScoreMatrix(alignment->matrix, alignment, vector_len, ofs_AVX);
        ofs_AVX.close();
    }

    free(alignment->matrix);
    free(alignment_lin->matrix);
LT_END_TEST(NW_ASM_AVX_test)

LT_BEGIN_TEST(TestNW, NW_ASM_AVX512_test)
    int vector_len = 32;
    int s_len1 = 32 + rand() % 100;
    int s_len2 = 32 + rand() % 100;
    

    char *s1 = random_seq(s_len1);
    char *s2 = random_seq(s_len2);
    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string(s1);
    alignment->sequence_2 = new_Sequence_from_string(s2);

    // Actualizo valor de s_len para considerar el -
    s_len1 = alignment->sequence_1->length;
    s_len2 = alignment->sequence_2->length;

    NW::C::with_logic_SSE(*alignment_lin,vector_len, true);
    
    NW_ASM_AVX512(alignment, true);

    bool score = (alignment->result->score == alignment_lin->result->score);
    LT_CHECK(score);

    bool valid_seqs = valid_alignment(*alignment_lin) && valid_alignment(*alignment);

    LT_CHECK(valid_seqs);

    bool matrix_ok = true;
    bool position_ok;
    for(int i = 0 ; i < s_len2 ; i++){
        for(int j = 0 ; j < s_len1 ; j++){
            position_ok = get_score_SSE(alignment->matrix,s_len1,i,j,vector_len) == get_score_SSE(alignment_lin->matrix,s_len1,i,j,vector_len);
            matrix_ok &= position_ok;
            if (!position_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_SSE(alignment_lin->matrix,s_len1,i,j,32));
                DBG(get_score_SSE(alignment->matrix,s_len1,i,j,32));   
            }
        }
    }
    
    LT_CHECK( matrix_ok );

    if (!matrix_ok || !valid_seqs || !score) {
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score AVX: %d\n", alignment->result->score);

        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment AVX512:\n");
        printf("%s\n", alignment->result->sequence_1->sequence);
        printf("%s\n", alignment->result->sequence_2->sequence);

        ofstream ofs_withLogicAVX512("NW_ASM_score_matrix_logic.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_lin->matrix, alignment_lin, vector_len, ofs_withLogicAVX512);
        ofs_withLogicAVX512.close();
        ofstream ofs_AVX512("NW_ASM_score_matrix_AVX512.txt", std::ofstream::trunc);
        printScoreMatrix(alignment->matrix, alignment, vector_len, ofs_AVX512);
        ofs_AVX512.close();
    }

    free(alignment->matrix);
    free(alignment_lin->matrix);
LT_END_TEST(NW_ASM_AVX512_test)

//===================================Smith waterman=====================================
LT_BEGIN_SUITE(TestSW)

void set_up() {
    printf("Ejecutando test para SW...\n");
}
void tear_down() {}

LT_END_SUITE(TestSW)

LT_BEGIN_TEST(TestSW, SW_C_x5)
    int vector_len = 32;
    int s_len1 = vector_len + rand() % 100;
    int s_len2 = vector_len + rand() % 100;
    char *s1 = random_seq(s_len1);
    char *s2 = random_seq(s_len2);
    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment_logic = new_alignment();
    alignment_logic->sequence_1 = new_Sequence_from_string(s1);
    alignment_logic->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment_sse = new_alignment();
    alignment_sse->sequence_1 = new_Sequence_from_string(s1);
    alignment_sse->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment_avx = new_alignment();
    alignment_avx->sequence_1 = new_Sequence_from_string(s1);
    alignment_avx->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment_avx512 = new_alignment();
    alignment_avx512->sequence_1 = new_Sequence_from_string(s1);
    alignment_avx512->sequence_2 = new_Sequence_from_string(s2);

    // Actualizo valor de s_len para considerar el -
    s_len1 = alignment_avx->sequence_1->length;
    s_len2 = alignment_avx->sequence_2->length;

    SW::C::LIN::SW(*alignment_lin, true);
    SW::C::with_logic_SSE(*alignment_logic, vector_len, true);
    SW::C::SSE::SW(*alignment_sse, true);
    SW::C::AVX::SW(*alignment_avx, true);
    SW::C::AVX512::SW(*alignment_avx512, true);

    bool score = (alignment_avx512->result->score == alignment_lin->result->score) &&
                 (alignment_avx->result->score == alignment_lin->result->score) &&
                 (alignment_sse->result->score == alignment_lin->result->score) &&
                 (alignment_logic->result->score == alignment_lin->result->score);

    LT_CHECK(score);

    bool valid_seqs = valid_alignment(*alignment_lin) &&
                      valid_alignment(*alignment_logic) &&
                      valid_alignment(*alignment_sse) &&
                      valid_alignment(*alignment_avx) &&
                      valid_alignment(*alignment_avx512);

    LT_CHECK(valid_seqs);

    bool matrix_ok = true;
    bool position_ok;
    for(int i = 0 ; i < s_len2 ; i++){
        for(int j = 0 ; j < s_len1 ; j++){
            position_ok = get_score_SSE(alignment_avx512->matrix,s_len1,i,j,32) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8) &&
                          get_score_SSE(alignment_avx->matrix,s_len1,i,j,16) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8) &&
                          get_score_SSE(alignment_sse->matrix,s_len1,i,j,8) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8) &&
                          get_score_SSE(alignment_logic->matrix,s_len1,i,j,32) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8);
            matrix_ok &= position_ok;
            if (!position_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_LIN(alignment_lin->matrix,s_len1,i,j,8));
                DBG(get_score_SSE(alignment_logic->matrix,s_len1,i,j,32));
                DBG(get_score_SSE(alignment_sse->matrix,s_len1,i,j,8));
                DBG(get_score_SSE(alignment_avx->matrix,s_len1,i,j,16));  
                DBG(get_score_SSE(alignment_avx512->matrix,s_len1,i,j,32));  
            }
        }
    }
    
    LT_CHECK( matrix_ok );
    
    if (!matrix_ok || !valid_seqs || !score) {
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score Logic: %d\n", alignment_logic->result->score);
        printf("Score SSE: %d\n", alignment_sse->result->score);
        printf("Score AVX: %d\n", alignment_avx->result->score);

        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment Logic:\n");
        printf("%s\n", alignment_logic->result->sequence_1->sequence);
        printf("%s\n", alignment_logic->result->sequence_2->sequence);

        printf("Alignment SSE:\n");
        printf("%s\n", alignment_sse->result->sequence_1->sequence);
        printf("%s\n", alignment_sse->result->sequence_2->sequence);
        
        printf("Alignment AVX:\n");
        printf("%s\n", alignment_avx->result->sequence_1->sequence);
        printf("%s\n", alignment_avx->result->sequence_2->sequence);

        printf("Alignment AVX512:\n");
        printf("%s\n", alignment_avx512->result->sequence_1->sequence);
        printf("%s\n", alignment_avx512->result->sequence_2->sequence);

        ofstream ofs_logic("SW_C_score_matrix_logic.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_logic->matrix, alignment_logic, 32, ofs_logic);
        ofs_logic.close();
        
        ofstream ofs_SSE("SW_C_score_matrix_SSE.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_sse->matrix, alignment_sse, 8, ofs_SSE);
        ofs_SSE.close();

        ofstream ofs_AVX("SW_C_score_matrix_AVX.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_avx->matrix, alignment_avx, 16, ofs_AVX);
        ofs_AVX.close();
        
        ofstream ofs_AVX512("SW_C_score_matrix_AVX512.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_avx512->matrix, alignment_avx512, 32, ofs_AVX512);
        ofs_AVX512.close();

    }

    free(alignment_avx512->matrix);
    free(alignment_avx->matrix);
    free(alignment_sse->matrix);
    free(alignment_logic->matrix);
    free(alignment_lin->matrix);
LT_END_TEST(SW_C_x5)

LT_BEGIN_TEST(TestSW, SW_ASM_SSE_test)
    int s_len1 = 8 + rand() % 100;
    int s_len2 = 8 + rand() % 100;
    char *s1 = random_seq(s_len1);
    char *s2 = random_seq(s_len2);

    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string(s1);
    alignment->sequence_2 = new_Sequence_from_string(s2);

    // Actualizo valor de s_len para considerar el -
    s_len1 = alignment->sequence_1->length;
    s_len2 = alignment->sequence_2->length;

    SW::C::LIN::SW(*alignment_lin, true);
    
    SW_ASM_SSE(alignment, true);

    bool score = (alignment->result->score == alignment_lin->result->score);
    LT_CHECK(score);

    bool valid_seqs = valid_alignment(*alignment_lin) && valid_alignment(*alignment);

    LT_CHECK(valid_seqs);

    bool matrix_ok = true;
    bool position_ok;
    for(int i = 0 ; i < s_len2 ; i++){
        for(int j = 0 ; j < s_len1 ; j++){
            position_ok = get_score_SSE(alignment->matrix,s_len1,i,j,8) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,8);
            matrix_ok &= position_ok;
            if (!position_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_SSE(alignment_lin->matrix,s_len1,i,j,8));
                DBG(get_score_SSE(alignment->matrix,s_len1,i,j,8));   
            }
        }
    }
    
    LT_CHECK( matrix_ok );

    if (!matrix_ok || !valid_seqs || !score) {
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score SSE: %d\n", alignment->result->score);

        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment SSE:\n");
        printf("%s\n", alignment->result->sequence_1->sequence);
        printf("%s\n", alignment->result->sequence_2->sequence);

        ofstream ofs_SSE("SW_ASM_score_matrix_SSE.txt", std::ofstream::trunc);
        printScoreMatrix(alignment->matrix, alignment, 8, ofs_SSE);
        ofs_SSE.close();
    }

    free(alignment->matrix);
    free(alignment_lin->matrix);
LT_END_TEST(SW_ASM_SSE_test)


LT_BEGIN_TEST(TestSW, SW_ASM_AVX_test)
    int s_len1 = 16 + rand() % 100;
    int s_len2 = 16 + rand() % 100;
    int vector_len = 16;
    char *s1 = random_seq(s_len1);
    char *s2 = random_seq(s_len2);

    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string(s1);
    alignment->sequence_2 = new_Sequence_from_string(s2);

    // Actualizo valor de s_len para considerar el -
    s_len1 = alignment->sequence_1->length;
    s_len2 = alignment->sequence_2->length;


    //SW::SW_C_withLogicSSE(*alignment_lin,16, true);
    SW::C::LIN::SW(*alignment_lin,true);
    
    SW_ASM_AVX(alignment, true);

    bool score = (alignment->result->score == alignment_lin->result->score);
    LT_CHECK(score);

    bool valid_seqs = valid_alignment(*alignment_lin) && valid_alignment(*alignment);

    LT_CHECK(valid_seqs);

    bool matrix_ok = true;
    bool position_ok;
    for(int i = 0 ; i < s_len2 ; i++){
        for(int j = 0 ; j < s_len1 ; j++){
            position_ok = get_score_SSE(alignment->matrix,s_len1,i,j,vector_len) == get_score_LIN(alignment_lin->matrix,s_len1,i,j,vector_len);
            matrix_ok &= position_ok;
             if (!position_ok){	
             printf("Score matrices differ at: \n");	
             DBG(i);	
             DBG(j);	
             DBG(get_score_LIN(alignment_lin->matrix,s_len1,i,j,8));	
             DBG(get_score_SSE(alignment->matrix,s_len1,i,j,8));   	
            }	
        }	
    }	
    	
    LT_CHECK( matrix_ok );	
    if (!matrix_ok || !valid_seqs || !score) {	
        printf("Score LIN: %d\n", alignment_lin->result->score);	
        printf("Score AVX: %d\n", alignment->result->score);
        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment AVX:\n");
        printf("%s\n", alignment->result->sequence_1->sequence);
        printf("%s\n", alignment->result->sequence_2->sequence);

        ofstream ofs_AVX("SW_ASM_score_matrix_AVX.txt", std::ofstream::trunc);
        printScoreMatrix(alignment->matrix, alignment, vector_len, ofs_AVX);
        ofs_AVX.close();
    }

    free(alignment->matrix);
    free(alignment_lin->matrix);
LT_END_TEST(SW_ASM_AVX_test)


LT_BEGIN_TEST(TestSW, SW_ASM_AVX512_test)
    int vector_len = 32;
    int s_len1 = 32 + rand() % 100;
    int s_len2 = 32 + rand() % 100;
    

    char *s1 = random_seq(s_len1);
    char *s2 = random_seq(s_len2);

    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string(s1);
    alignment->sequence_2 = new_Sequence_from_string(s2);

    // Actualizo valor de s_len para considerar el -
    s_len1 = alignment->sequence_1->length;
    s_len2 = alignment->sequence_2->length;

    SW::C::with_logic_SSE(*alignment_lin,vector_len, true);
    
    SW_ASM_AVX512(alignment, true);

    bool score = (alignment->result->score == alignment_lin->result->score);
    LT_CHECK(score);

    bool valid_seqs = valid_alignment(*alignment_lin) && valid_alignment(*alignment);

    LT_CHECK(valid_seqs);

    bool matrix_ok = true;
    bool position_ok;
    for(int i = 0 ; i < s_len2 ; i++){
        for(int j = 0 ; j < s_len1 ; j++){
            position_ok = get_score_SSE(alignment->matrix,s_len1,i,j,vector_len) == get_score_SSE(alignment_lin->matrix,s_len1,i,j,vector_len);
            matrix_ok &= position_ok;
            if (!position_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_SSE(alignment_lin->matrix,s_len1,i,j,32));
                DBG(get_score_SSE(alignment->matrix,s_len1,i,j,32));   
            }
        }
    }
    
    LT_CHECK( matrix_ok );

    if (!matrix_ok || !valid_seqs || !score) {
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score AVX: %d\n", alignment->result->score);

        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment AVX512:\n");
        printf("%s\n", alignment->result->sequence_1->sequence);
        printf("%s\n", alignment->result->sequence_2->sequence);

        ofstream ofs_withLogicAVX512("SW_ASM_score_matrix_logic.txt", std::ofstream::trunc);
        printScoreMatrix(alignment_lin->matrix, alignment_lin, vector_len, ofs_withLogicAVX512);
        ofs_withLogicAVX512.close();
        ofstream ofs_AVX512("SW_ASM_score_matrix_AVX512.txt", std::ofstream::trunc);
        printScoreMatrix(alignment->matrix, alignment, vector_len, ofs_AVX512);
        ofs_AVX512.close();
    }

    free(alignment->matrix);
    free(alignment_lin->matrix);
LT_END_TEST(SW_ASM_AVX512_test)


// Ejecutar tests
LT_BEGIN_AUTO_TEST_ENV()
    // Set random seed
    srand( (unsigned) time(NULL) * getpid());
    // Run
    AUTORUN_TESTS()
LT_END_AUTO_TEST_ENV()
