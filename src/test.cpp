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

//===========================Needleman wunsch===========================================
LT_BEGIN_SUITE(TestNW)

void set_up() {
    printf("Ejecutando test para NW...\n");
}
void tear_down() {}

LT_END_SUITE(TestNW)

LT_BEGIN_TEST(TestNW, NW_C_LIN)

    Alignment* alignment = new_alignment();
    alignment->parameters->gap = -1;
    alignment->parameters->match = 1;
    alignment->parameters->missmatch = -1;
    alignment->sequence_1 = new_Sequence_from_string((char*) "GCATGCU");
    alignment->sequence_2 = new_Sequence_from_string((char*) "GATTACA");
    
    NW::NW_C_LIN(*alignment, true);

    bool res1 = strcmp(alignment->result->sequence_2->sequence, "-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCAT-GCU") == 0;
    bool res2 = strcmp(alignment->result->sequence_2->sequence,"-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCA-TGCU") == 0;
    bool res3 = strcmp(alignment->result->sequence_2->sequence,"-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCATG-CU") == 0;

    LT_CHECK((res1 || res2 || res3) && alignment->result->score == 0);

    short matrix_nw[8][8] = {{0,-1,-2,-3,-4,-5,-6,-7},
                              {-1,1,0,-1,-2,-3,-4,-5},
                              {-2,0,0,1,0,-1,-2,-3},
                              {-3,-1,-1,0,2,1,0,-1},
                              {-4,-2,-2,-1,1,1,0,-1},
                              {-5,-3,-3,-1,0,0,0,-1},
                              {-6,-4,-2,-2,-1,-1,1,0},
                              {-7,-5,-3,-1,-2,-2,0,0}};

    
    bool check_nw = check_scr_matrix_manual((short **)matrix_nw, alignment, get_score_LIN);
    
    LT_CHECK( check_nw );
    free(alignment->matrix);
LT_END_TEST(NW_C_LIN)

LT_BEGIN_TEST(TestNW, NW_C_LIN2)

    Alignment* alignment = new_alignment();
    alignment->parameters->gap = -2;
    alignment->parameters->match = 1;
    alignment->parameters->missmatch = -1;
    alignment->sequence_1 = new_Sequence_from_string((char*) "TGGTG");
    alignment->sequence_2 = new_Sequence_from_string((char*) "ATCGT");
    
    NW::NW_C_LIN(*alignment, true);

    bool res = strcmp(alignment->result->sequence_2->sequence, "-ATCGT-") == 0 && strcmp(alignment->result->sequence_1->sequence, "--TGGTG") == 0;

    LT_CHECK(res && alignment->result->score == -2);

    short matrix_nw[6][6] = {{0,-2,-4,-6,-8,-10},
                             {-2,-1,-3,-5,-7,-9},
                             {-4,-1,-2,-4,-4,-6},
                             {-6,-3,-2,-3,-5,-5},
                             {-8,-5,-2,-1,-3,-4},
                             {-10,-7,-4,-3,0,-2}};
    
    bool check_nw = check_scr_matrix_manual((short **)matrix_nw, alignment, get_score_LIN);
    
    LT_CHECK( check_nw );
    free(alignment->matrix);
LT_END_TEST(NW_C_LIN2)

// LT_BEGIN_TEST(TestNW, NW_C_withLogicSSE)
//     Alignment* alignment = new_alignment();
//     alignment->sequence_1 = new_Sequence_from_string((char*) "GCATGCU");
//     alignment->sequence_2 = new_Sequence_from_string((char*) "GATTACA");
    
//     NW::NW_C_withLogicSSE(*alignment, true);

//     bool res1 = strcmp(alignment->result->sequence_2->sequence, "-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCAT-GCU") == 0;
//     bool res2 = strcmp(alignment->result->sequence_2->sequence,"-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCA-TGCU") == 0;
//     bool res3 = strcmp(alignment->result->sequence_2->sequence,"-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCATG-CU") == 0;

//     LT_CHECK((res1 || res2 || res3) && alignment->result->score == 0);

//     short matrix_nw[8][8] = {{0,-1,-2,-3,-4,-5,-6,-7},
//                               {-1,1,0,-1,-2,-3,-4,-5},
//                               {-2,0,0,1,0,-1,-2,-3},
//                               {-3,-1,-1,0,2,1,0,-1},
//                               {-4,-2,-2,-1,1,1,0,-1},
//                               {-5,-3,-3,-1,0,0,0,-1},
//                               {-6,-4,-2,-2,-1,-1,1,0},
//                               {-7,-5,-3,-1,-2,-2,0,0}};
    
    
//     bool check_nw = check_scr_matrix_manual((short **)matrix_nw, alignment, get_score_SSE);
    
//     LT_CHECK( check_nw );
    
// LT_END_TEST(NW_C_withLogicSSE)

// LT_BEGIN_TEST(TestNW, NW_C_SSE)
//     Alignment* alignment = new_alignment();
//     alignment->sequence_1 = new_Sequence_from_string((char*) "GCATGCU");
//     alignment->sequence_2 = new_Sequence_from_string((char*) "GATTACA");
    
//     NW::NW_C_SSE(*alignment, true);

//     bool res1 = strcmp(alignment->result->sequence_2->sequence, "-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCAT-GCU") == 0;
//     bool res2 = strcmp(alignment->result->sequence_2->sequence,"-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCA-TGCU") == 0;
//     bool res3 = strcmp(alignment->result->sequence_2->sequence,"-G-ATTACA") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GCATG-CU") == 0;

//     LT_CHECK((res1 || res2 || res3) && alignment->result->score == 0);

//     short matrix_res[8][8] = {{0,-1,-2,-3,-4,-5,-6,-7},
//                               {-1,1,0,-1,-2,-3,-4,-5},
//                               {-2,0,0,1,0,-1,-2,-3},
//                               {-3,-1,-1,0,2,1,0,-1},
//                               {-4,-2,-2,-1,1,1,0,-1},
//                               {-5,-3,-3,-1,0,0,0,-1},
//                               {-6,-4,-2,-2,-1,-1,1,0},
//                               {-7,-5,-3,-1,-2,-2,0,0}};

//     bool check_nw = check_scr_matrix_manual((short **)matrix_res, alignment, get_score_SSE);
    
//     LT_CHECK( check_nw );
    
// LT_END_TEST(NW_C_SSE)

//===================================Smith waterman=====================================
LT_BEGIN_SUITE(TestSW)

void set_up() {
    printf("Ejecutando test para SW...\n");
}
void tear_down() {}

LT_END_SUITE(TestSW)

LT_BEGIN_TEST(TestSW, SW_C_LIN)
    Alignment* alignment = new_alignment();
    alignment->parameters->gap = -2;
    alignment->parameters->match = 3;
    alignment->parameters->missmatch = -3;
    alignment->sequence_1 = new_Sequence_from_string((char*) "TGTTACGG");
    alignment->sequence_2 = new_Sequence_from_string((char*) "GGTTGACTA");
    
    SW::SW_C_LIN(*alignment, true);

    bool res = strcmp(alignment->result->sequence_1->sequence, "-GTT-AC") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GTTGAC") == 0;
    
    LT_CHECK(res && alignment->result->score == 13);

    short matrix_sw[10][9] = {{0,0,0,0,0,0,0,0,0},
                              {0,0,3,1,0,0,0,3,3},
                              {0,0,3,1,0,0,0,3,6},
                              {0,3,1,6,4,2,0,1,4},
                              {0,3,1,4,9,7,5,3,2},
                              {0,1,6,4,7,6,4,8,6},
                              {0,0,4,3,5,10,8,6,5},
                              {0,0,2,1,3,8,13,11,9},
                              {0,3,1,5,4,6,11,10,8},
                              {0,1,0,3,2,7,9,8,7}};
    
    bool check_sw = check_scr_matrix_manual((short **)matrix_sw, alignment, get_score_LIN);
    
    LT_CHECK( check_sw );

    free(alignment->matrix);   
LT_END_TEST(SW_C_LIN) 

LT_BEGIN_TEST(TestNW, SW_C_LIN2)

    Alignment* alignment = new_alignment();
    alignment->parameters->gap = -2;
    alignment->parameters->match = 1;
    alignment->parameters->missmatch = -1;
    alignment->sequence_1 = new_Sequence_from_string((char*) "TGGTG");
    alignment->sequence_2 = new_Sequence_from_string((char*) "ATCGT");
    
    SW::SW_C_LIN(*alignment, true);

    bool res = strcmp(alignment->result->sequence_2->sequence, "-GT") == 0 && strcmp(alignment->result->sequence_1->sequence, "-GT") == 0;

    LT_CHECK(res && alignment->result->score == 2);

    short matrix_sw[6][6] = {{0,0,0,0,0,0},
                             {0,0,0,0,0,0},
                             {0,1,0,0,1,0},
                             {0,0,0,0,0,0},
                             {0,0,1,1,0,1},
                             {0,1,0,0,2,0}};
    
    bool check_sw = check_scr_matrix_manual((short **)matrix_sw, alignment, get_score_LIN);

    LT_CHECK( check_sw );
    free(alignment->matrix);
LT_END_TEST(SW_C_LIN2)

// LT_BEGIN_TEST(TestSW, SW_C_withLogicSSE)
//      Alignment* alignment = new_alignment();
//     alignment->parameters->gap = -2;
//     alignment->parameters->match = 3;
//     alignment->parameters->missmatch = -3;
//     alignment->sequence_1 = new_Sequence_from_string((char*) "TGTTACGG");
//     alignment->sequence_2 = new_Sequence_from_string((char*) "GGTTGACTA");
    
//     SW::SW_C_withLogicSSE(*alignment, true);

//     bool res = strcmp(alignment->result->sequence_1->sequence, "-GTT-AC") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GTTGAC") == 0;
    
//     LT_CHECK(res && alignment->result->score == 13);

//     short matrix_sw[10][9] = {{0,0,0,0,0,0,0,0,0},
//                               {0,0,3,1,0,0,0,3,3},
//                               {0,0,3,1,0,0,0,3,6},
//                               {0,3,1,6,4,2,0,1,4},
//                               {0,3,1,4,9,7,5,3,2},
//                               {0,1,6,4,7,6,4,8,6},
//                               {0,0,4,3,5,10,8,6,5},
//                               {0,0,2,1,3,8,13,11,9},
//                               {0,3,1,5,4,6,11,10,8},
//                               {0,1,0,3,2,7,9,8,7}};

//     bool check_sw = check_scr_matrix_manual((short **)matrix_sw, alignment, get_score_SSE);
    
//     LT_CHECK( check_sw );

//     free(alignment->matrix);   
// LT_END_TEST(SW_C_withLogicSSE)

// LT_BEGIN_TEST(TestSW, SW_C_SSE)
//      Alignment* alignment = new_alignment();
//     alignment->parameters->gap = -2;
//     alignment->parameters->match = 3;
//     alignment->parameters->missmatch = -3;
//     alignment->sequence_1 = new_Sequence_from_string((char*) "TGTTACGG");
//     alignment->sequence_2 = new_Sequence_from_string((char*) "GGTTGACTA");
    
//     SW::SW_C_SSE(*alignment, true);

//     bool res = strcmp(alignment->result->sequence_1->sequence, "-GTT-AC") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GTTGAC") == 0;
    
//     LT_CHECK(res && alignment->result->score == 13);

//     short matrix_sw[10][9] = {{0,0,0,0,0,0,0,0,0},
//                               {0,0,3,1,0,0,0,3,3},
//                               {0,0,3,1,0,0,0,3,6},
//                               {0,3,1,6,4,2,0,1,4},
//                               {0,3,1,4,9,7,5,3,2},
//                               {0,1,6,4,7,6,4,8,6},
//                               {0,0,4,3,5,10,8,6,5},
//                               {0,0,2,1,3,8,13,11,9},
//                               {0,3,1,5,4,6,11,10,8},
//                               {0,1,0,3,2,7,9,8,7}};
    
//     bool check_sw = check_scr_matrix_manual((short **)matrix_sw, alignment, get_score_SSE);
    
//     LT_CHECK( check_sw );

//     free(alignment->matrix);   
// LT_END_TEST(SW_C_SSE)

// LT_BEGIN_TEST(TestSW, SW_C_SSE_EQUAL_STRINGS)
//      Alignment* alignment = new_alignment();
//     alignment->parameters->gap = -2;
//     alignment->parameters->match = 3;
//     alignment->parameters->missmatch = -3;
//     alignment->sequence_1 = new_Sequence_from_string((char*) "TTTTTTTTT");
//     alignment->sequence_2 = new_Sequence_from_string((char*) "TTTTTTTTT");
    
//     SW::SW_C_SSE(*alignment, true);

//     bool res = strcmp(alignment->result->sequence_1->sequence, "-TTTTTTTTT") == 0 && strcmp(alignment->result->sequence_2->sequence, "-TTTTTTTTT") == 0;
//     printf("%s\n",alignment->result->sequence_1->sequence);
//     printf("%s\n",alignment->result->sequence_2->sequence);
//     LT_CHECK(res);

//     free(alignment->matrix);   
// LT_END_TEST(SW_C_SSE_EQUAL_STRINGS)


// LT_BEGIN_TEST(TestSW, SW_C_SSE_TEST1)
//      Alignment* alignment = new_alignment();
//     alignment->parameters->gap = -2;
//     alignment->parameters->match = 3;
//     alignment->parameters->missmatch = -3;
//     alignment->sequence_1 = new_Sequence_from_string((char*) "TTTTTATT");
//     alignment->sequence_2 = new_Sequence_from_string((char*) "TTTTTTTATT");
    
//     SW::SW_C_SSE(*alignment, true);

//     bool res = strcmp(alignment->result->sequence_1->sequence, "-TTTTTATT") == 0 && strcmp(alignment->result->sequence_2->sequence, "-TTTTTATT") == 0;
//     printf("%s\n",alignment->result->sequence_1->sequence);
//     printf("%s\n",alignment->result->sequence_2->sequence);
//     LT_CHECK(res);

//     free(alignment->matrix);   
// LT_END_TEST(SW_C_SSE_TEST1)

// LT_BEGIN_TEST(TestSW, SW_C_SSE_FINALTEST)
//      Alignment* alignment = new_alignment();
//     alignment->parameters->gap = -2;
//     alignment->parameters->match = 3;
//     alignment->parameters->missmatch = -3;
//     alignment->sequence_1 = new_Sequence_from_string((char*) "TGTTACGG");
//     alignment->sequence_2 = new_Sequence_from_string((char*) "GGTTGACTA");
    
//     SW::SW_C_SSE(*alignment, true);

//     bool res = strcmp(alignment->result->sequence_1->sequence, "-GTT-AC") == 0 && strcmp(alignment->result->sequence_2->sequence, "-GTTGAC") == 0;
//     printf("%s\n",alignment->result->sequence_1->sequence);
//     printf("%s\n",alignment->result->sequence_2->sequence);
    
//     LT_CHECK( res );

//     free(alignment->matrix);   
// LT_END_TEST(SW_C_SSE_FINALTEST)


// Ejecutar tests
LT_BEGIN_AUTO_TEST_ENV()
    AUTORUN_TESTS()
LT_END_AUTO_TEST_ENV()
