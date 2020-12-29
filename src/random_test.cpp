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


#define DBG(x) cerr << #x << " = " << (x) <<"\n"
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

LT_BEGIN_TEST(TestNW, NW_C_SSE)
    int s_len = 5 + rand() % 96;
    char *s1 = random_seq(s_len);
    char *s2 = random_seq(s_len);

    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string(s1);
    alignment->sequence_2 = new_Sequence_from_string(s2);

    NW::NW_C_LIN(*alignment_lin, true);
    
    NW::NW_C_SSE(*alignment, true);

    bool score = alignment->result->score == alignment_lin->result->score;

    LT_CHECK(score);

    if (!score){
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score SSE: %d\n", alignment->result->score);
    }

    bool seqs = strcmp(alignment->result->sequence_1->sequence, alignment_lin->result->sequence_1->sequence) == 0 &&
                strcmp(alignment->result->sequence_2->sequence, alignment_lin->result->sequence_2->sequence) == 0;

    LT_CHECK(seqs);

    if (!seqs) {
        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment SSE:\n");
        printf("%s\n", alignment->result->sequence_1->sequence);
        printf("%s\n", alignment->result->sequence_2->sequence);
    }

    bool matrix_ok = true;
    for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
            matrix_ok &= get_score_SSE(alignment->matrix->matrix,s_len,i,j,8) == get_score_LIN(alignment_lin->matrix->matrix,s_len,i,j,8);
            if (!matrix_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_SSE(alignment->matrix->matrix,s_len,i,j,8));
                DBG(get_score_LIN(alignment_lin->matrix->matrix,s_len,i,j,8));
                break;       
            }
        }
        if (!matrix_ok)
            break; 
    }
    
    LT_CHECK( matrix_ok );

    free(alignment->matrix->matrix);
    free(alignment_lin->matrix->matrix);
LT_END_TEST(NW_C_SSE)

//===================================Smith waterman=====================================
LT_BEGIN_SUITE(TestSW)

void set_up() {
    printf("Ejecutando test para SW...\n");
}
void tear_down() {}

LT_END_SUITE(TestSW)

LT_BEGIN_TEST(TestSW, SW_C_SSE)
    int s_len = 5 + rand() % 96;
    char *s1 = random_seq(s_len);
    char *s2 = random_seq(s_len);

    printf("Input sequence 1: %s\n", s1);
    printf("Input sequence 2: %s\n", s2);

    Alignment* alignment_lin = new_alignment();
    alignment_lin->sequence_1 = new_Sequence_from_string(s1);
    alignment_lin->sequence_2 = new_Sequence_from_string(s2);

    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string(s1);
    alignment->sequence_2 = new_Sequence_from_string(s2);

    SW::SW_C_LIN(*alignment_lin, true);
    
    SW::SW_C_SSE(*alignment, true);

    bool score = alignment->result->score == alignment_lin->result->score;

    LT_CHECK(score);

    if (!score){
        printf("Score LIN: %d\n", alignment_lin->result->score);
        printf("Score SSE: %d\n", alignment->result->score);
    }

    bool seqs = strcmp(alignment->result->sequence_1->sequence, alignment_lin->result->sequence_1->sequence) == 0 &&
                strcmp(alignment->result->sequence_2->sequence, alignment_lin->result->sequence_2->sequence) == 0;
    
    LT_CHECK(seqs);

    if (!seqs) {
        printf("Alignment LIN:\n");
        printf("%s\n", alignment_lin->result->sequence_1->sequence);
        printf("%s\n", alignment_lin->result->sequence_2->sequence);

        printf("Alignment SSE:\n");
        printf("%s\n", alignment->result->sequence_1->sequence);
        printf("%s\n", alignment->result->sequence_2->sequence);
    }

    bool matrix_ok = true;
    for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
            matrix_ok &= get_score_SSE(alignment->matrix->matrix,s_len,i,j,8) == get_score_LIN(alignment_lin->matrix->matrix,s_len,i,j,8);
            if (!matrix_ok){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(get_score_SSE(alignment->matrix->matrix,s_len,i,j,8));
                DBG(get_score_LIN(alignment_lin->matrix->matrix,s_len,i,j,8));
                break;       
            }
        }
        if (!matrix_ok)
            break; 
    }
    
    LT_CHECK( matrix_ok );

    free(alignment->matrix->matrix);
    free(alignment_lin->matrix->matrix);
LT_END_TEST(SW_C_SSE)

// Ejecutar tests
LT_BEGIN_AUTO_TEST_ENV()
    // Set random seed
    srand( (unsigned) time(NULL) * getpid());
    // Run
    AUTORUN_TESTS()
LT_END_AUTO_TEST_ENV()
