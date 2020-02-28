CXX=g++
CXXFLAGS = -ggdb -Wall -Wno-unused-parameter -Wextra -no-pie -pedantic -m64 -O0 -march=native

ASM = nasm
ASMFLAGS = -felf64 -g -F dwarf

BUILD_DIR = ../build

ALGORITMOS = NW SW

ALGORITMOS_OBJ = $(addsuffix _C.o, $(ALGORITMOS)) $(addsuffix _AVX.o, $(ALGORITMOS)) $(addsuffix _SSE.o, $(ALGORITMOS)) $(addsuffix _LIN.o, $(ALGORITMOS))
ALGORITMOS_OBJ_CON_PATH = $(addprefix  $(BUILD_DIR)/, $(ALGORITMOS_OBJ))

all: cli

cli: $(BUILD_DIR)/cli.o $(ALGORITMOS_OBJ_CON_PATH) $(BUILD_DIR)/Types.o
	$(CXX) $(CXXFLAGS) -o ../$@ $(BUILD_DIR)/cli.o $(ALGORITMOS_OBJ_CON_PATH) $(BUILD_DIR)/Types.o

$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(BUILD_DIR)/%.o: %.asm
	mkdir -p $(BUILD_DIR)
	$(ASM) $(ASMFLAGS) -o $@ $<

clean:
	rm -f $(BUILD_DIR)/*
	rm -f ../cli