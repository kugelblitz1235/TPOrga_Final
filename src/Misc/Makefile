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

all: cli

cli: $(BUILD_DIR)/cli.o $(ALGORITMOS_OBJ_CON_PATH) $(AUX_OBJ_CON_PATH)
	$(CXX) $(CXXFLAGS) -o ../$@ $(BUILD_DIR)/cli.o $(ALGORITMOS_OBJ_CON_PATH) $(AUX_OBJ_CON_PATH)

$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(BUILD_DIR)/%.o: %.asm
	mkdir -p $(BUILD_DIR)
	$(ASM) $(ASMFLAGS) -o $@ $<

# $(BUILD_DIR)/Utility.o: Utility.cpp
# 	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -f $(BUILD_DIR)/*
	rm -f ../cli