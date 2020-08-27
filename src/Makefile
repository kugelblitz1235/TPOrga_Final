CXX=g++
CXXFLAGS = -ggdb -Wall -Wno-unused-parameter -Wextra -no-pie -pedantic -m64 -O0 -march=native

ASM = nasm
ASMFLAGS = -felf64 -g -F dwarf

BUILD_DIR = ../build

MISC_FILES =  Types.o  Utility.o FASTA.o
MISC_FILES_CON_PATH = $(addprefix  $(BUILD_DIR)/, $(MISC_FILES))
NEEDLMAN_C_FILES = NW_C.o
NEEDLMAN_C_FILES_CON_PATH = $(addprefix  $(BUILD_DIR)/, $(NEEDLMAN_C_FILES))
SMITH_C_FILES = SW_C.o
SMITH_C_FILES_CON_PATH = $(addprefix  $(BUILD_DIR)/, $(SMITH_C_FILES))
NEEDLMAN_ASM_FILES = NW_AVX.o  NW_LIN.o  NW_SSE.o
NEEDLMAN_ASM_FILES_CON_PATH = $(addprefix  $(BUILD_DIR)/, $(NEEDLMAN_ASM_FILES))
SMITH_ASM_FILES = SW_AVX.o  SW_LIN.o  SW_SSE.o 
SMITH_ASM_FILES_CON_PATH = $(addprefix  $(BUILD_DIR)/, $(SMITH_ASM_FILES))



all: cli

cli: $(BUILD_DIR)/cli.o $(MISC_FILES_CON_PATH) $(NEEDLMAN_C_FILES_CON_PATH) $(SMITH_C_FILES_CON_PATH) $(NEEDLMAN_ASM_FILES_CON_PATH) $(SMITH_ASM_FILES_CON_PATH)
	$(CXX) $(CXXFLAGS) -o ../$@ $(BUILD_DIR)/cli.o $(MISC_FILES_CON_PATH) $(NEEDLMAN_C_FILES_CON_PATH) $(SMITH_C_FILES_CON_PATH) $(NEEDLMAN_ASM_FILES_CON_PATH) $(SMITH_ASM_FILES_CON_PATH)

$(BUILD_DIR)/cli.o: cli.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<
	
$(filter $(BUILD_DIR)/%.o, $(MISC_FILES_CON_PATH)): $(BUILD_DIR)/%.o: Misc/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<
	
$(filter $(BUILD_DIR)/%.o, $(NEEDLMAN_C_FILES_CON_PATH)): $(BUILD_DIR)/%.o: Needleman_Wunsch/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<
	
$(filter $(BUILD_DIR)/%.o, $(SMITH_C_FILES_CON_PATH)): $(BUILD_DIR)/%.o: Smith_Waterman/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(filter $(BUILD_DIR)/%.o, $(NEEDLMAN_ASM_FILES_CON_PATH)): $(BUILD_DIR)/%.o: Needleman_Wunsch/%.asm
	mkdir -p $(BUILD_DIR)
	$(ASM) $(ASMFLAGS) -o $@ $<

$(filter $(BUILD_DIR)/%.o, $(SMITH_ASM_FILES_CON_PATH)): $(BUILD_DIR)/%.o: Smith_Waterman/%.asm
	mkdir -p $(BUILD_DIR)
	$(ASM) $(ASMFLAGS) -o $@ $<

clean:
	rm -f $(BUILD_DIR)/*
	rm -f ../cli

test: $(BUILD_DIR)/test.o $(MISC_FILES_CON_PATH) $(NEEDLMAN_C_FILES_CON_PATH) $(SMITH_C_FILES_CON_PATH) $(NEEDLMAN_ASM_FILES_CON_PATH) $(SMITH_ASM_FILES_CON_PATH)
	$(CXX) $(CXXFLAGS) -o ../$@ $(BUILD_DIR)/test.o $(MISC_FILES_CON_PATH) $(NEEDLMAN_C_FILES_CON_PATH) $(SMITH_C_FILES_CON_PATH) $(NEEDLMAN_ASM_FILES_CON_PATH) $(SMITH_ASM_FILES_CON_PATH)
	../test

$(BUILD_DIR)/test.o: test.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<