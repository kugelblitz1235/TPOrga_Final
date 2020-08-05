CXX=g++
CXXFLAGS = -ggdb -Wall -Wno-unused-parameter -Wextra -no-pie -pedantic -m64 -O0 -march=native

ASM = nasm
ASMFLAGS = -felf64 -g -F dwarf

BUILD_DIR = ../build

ALGORITMOS = NW SW

ALGORITMOS_OBJ = $(addsuffix _C.o, $(ALGORITMOS)) $(addsuffix _AVX.o, $(ALGORITMOS)) $(addsuffix _SSE.o, $(ALGORITMOS)) $(addsuffix _LIN.o, $(ALGORITMOS))
ALGORITMOS_OBJ_CON_PATH = $(addprefix  $(BUILD_DIR)/, $(ALGORITMOS_OBJ))

AUX = Types Utility
AUX_OBJ = $(addsuffix .o, $(AUX))
AUX_OBJ_CON_PATH = $(addprefix  $(BUILD_DIR)/, $(AUX_OBJ))

MISC_FILES =  Types.o  Utility.o
NEEDLMAN_FILES = NW_AVX.o  NW_C.o  NW_LIN.o  NW_SSE.o
SMITH_FILES = SW_AVX.o  SW_C.o  SW_LIN.o  SW_SSE.o 

all: cli

cli: $(BUILD_DIR)/cli.o $(ALGORITMOS_OBJ_CON_PATH) $(AUX_OBJ_CON_PATH)
	$(CXX) $(CXXFLAGS) -o ../$@ $(BUILD_DIR)/cli.o $(ALGORITMOS_OBJ_CON_PATH) $(AUX_OBJ_CON_PATH)

$(BUILD_DIR)/cli.o: cli.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(filter $(BUILD_DIR)/%.o, $(MISC_FILES)): Misc/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c Misc/$<

$(filter $(BUILD_DIR)/%.o, $(NEEDLMAN_FILES)): Needleman_Wunsch/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c Needleman_Wunsch/$<

$(filter $(BUILD_DIR)/%.o, $(SMITH_FILES)): Smith_Waterman/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c Smith_Waterman/$<

$(filter $(BUILD_DIR)/%.o, $(NEEDLMAN_FILES)): Needleman_Wunsch/%.asm
	mkdir -p $(BUILD_DIR)
	$(ASM) $(ASMFLAGS) -o $@ Needleman_Wunsch/$<

$(filter $(BUILD_DIR)/%.o, $(SMITH_FILES)): Smith_Waterman/%.asm
	mkdir -p $(BUILD_DIR)
	$(ASM) $(ASMFLAGS) -o $@ Smith_Waterman/$<


# $(BUILD_DIR)/Utility.o: Utility.cpp
# 	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -f $(BUILD_DIR)/*
	rm -f ../cli