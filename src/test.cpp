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

LT_BEGIN_SUITE(TestNW)

void set_up() {
    printf("Ejecutando tests...\n");
}
void tear_down() {}

LT_END_SUITE(TestNW)

LT_BEGIN_TEST(TestNW, C_LIN)
    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string((char*) "GCATGCU");
    alignment->sequence_2 = new_Sequence_from_string((char*) "GATTACA");
    
    NW_C_LIN(*alignment, true);

    bool res1 = strcmp(alignment->result->sequence_1->sequence, "-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCAT_GCU") == 0;
    bool res2 = strcmp(alignment->result->sequence_1->sequence,"-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCA_TGCU") == 0;
    bool res3 = strcmp(alignment->result->sequence_1->sequence,"-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCATG_CU") == 0;

    LT_CHECK((res1 || res2 || res3) && alignment->result->score == 0);

    short matrix_res[8][8] = {{0,-1,-2,-3,-4,-5,-6,-7},
                              {-1,1,0,-1,-2,-3,-4,-5},
                              {-2,0,0,1,0,-1,-2,-3},
                              {-3,-1,-1,0,2,1,0,-1},
                              {-4,-2,-2,-1,1,1,0,-1},
                              {-5,-3,-3,-1,0,0,0,-1},
                              {-6,-4,-2,-2,-1,-1,1,0},
                              {-7,-5,-3,-1,-2,-2,0,0}};
    bool equals = true;
    for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
            equals &= matrix_res[i][j] == get_score_LIN(alignment->matrix->matrix,8,i,j);
            if (!equals){
                DBG(i);
                DBG(j);
                DBG(matrix_res[i][j]);
                DBG(get_score_LIN(alignment->matrix->matrix,8,i,j));
                break;       
            }
        

        }
        if (!equals)
            break; 
    }
    
    LT_CHECK( equals );
    
LT_END_TEST(C_LIN)


LT_BEGIN_TEST(TestNW, C_withLogicSSE)
    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string((char*) "GCATGCU");
    alignment->sequence_2 = new_Sequence_from_string((char*) "GATTACA");
    
    NW_C_withLogicSSE(*alignment, true);

    bool res1 = strcmp(alignment->result->sequence_1->sequence, "-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCAT_GCU") == 0;
    bool res2 = strcmp(alignment->result->sequence_1->sequence,"-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCA_TGCU") == 0;
    bool res3 = strcmp(alignment->result->sequence_1->sequence,"-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCATG_CU") == 0;

    LT_CHECK((res1 || res2 || res3) && alignment->result->score == 0);

    short matrix_res[8][8] = {{0,-1,-2,-3,-4,-5,-6,-7},
                              {-1,1,0,-1,-2,-3,-4,-5},
                              {-2,0,0,1,0,-1,-2,-3},
                              {-3,-1,-1,0,2,1,0,-1},
                              {-4,-2,-2,-1,1,1,0,-1},
                              {-5,-3,-3,-1,0,0,0,-1},
                              {-6,-4,-2,-2,-1,-1,1,0},
                              {-7,-5,-3,-1,-2,-2,0,0}};
    bool equals = true;
    for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
            equals &= matrix_res[i][j] == get_score_SSE(alignment->matrix->matrix,8,i,j,8);
            if (!equals){
                DBG(i);
                DBG(j);
                DBG(matrix_res[i][j]);
                DBG(get_score_SSE(alignment->matrix->matrix,8,i,j,8));
                break;       
            }
        

        }
        if (!equals)
            break; 
    }
    
    LT_CHECK( equals );
    
LT_END_TEST(C_withLogicSSE)

LT_BEGIN_TEST(TestNW, C_SSE)
    Alignment* alignment = new_alignment();
    alignment->sequence_1 = new_Sequence_from_string((char*) "GCATGCU");
    alignment->sequence_2 = new_Sequence_from_string((char*) "GATTACA");
    
    NW_C_SSE(*alignment, true);

    bool res1 = strcmp(alignment->result->sequence_1->sequence, "-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCAT_GCU") == 0;
    bool res2 = strcmp(alignment->result->sequence_1->sequence,"-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCA_TGCU") == 0;
    bool res3 = strcmp(alignment->result->sequence_1->sequence,"-G_ATTACA") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GCATG_CU") == 0;

    LT_CHECK((res1 || res2 || res3) && alignment->result->score == 0);

    short matrix_res[8][8] = {{0,-1,-2,-3,-4,-5,-6,-7},
                              {-1,1,0,-1,-2,-3,-4,-5},
                              {-2,0,0,1,0,-1,-2,-3},
                              {-3,-1,-1,0,2,1,0,-1},
                              {-4,-2,-2,-1,1,1,0,-1},
                              {-5,-3,-3,-1,0,0,0,-1},
                              {-6,-4,-2,-2,-1,-1,1,0},
                              {-7,-5,-3,-1,-2,-2,0,0}};
    bool equals = true;
    for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
            equals &= matrix_res[i][j] == get_score_SSE(alignment->matrix->matrix,8,i,j,8);
            if (!equals){
                DBG(i);
                DBG(j);
                DBG(matrix_res[i][j]);
                DBG(get_score_SSE(alignment->matrix->matrix,8,i,j,8));
                break;       
            }
        

        }
        if (!equals)
            break; 
    }

    
    LT_CHECK( equals );
    
LT_END_TEST(C_SSE)

// Ejecutar tests
LT_BEGIN_AUTO_TEST_ENV()
    AUTORUN_TESTS()
LT_END_AUTO_TEST_ENV()