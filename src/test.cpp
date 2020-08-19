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
    
LT_END_TEST(C_LIN)


// Ejecutar tests
LT_BEGIN_AUTO_TEST_ENV()
    AUTORUN_TESTS()
LT_END_AUTO_TEST_ENV()