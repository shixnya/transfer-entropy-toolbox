BIN_DIR = bin

all: te_block te_block_fixed te_block_1 example

te_block_1: te_block_1.cpp
	mkdir -p $(BIN_DIR)
	g++ -O2 -Wall -lboost_program_options -lboost_thread -o $(BIN_DIR)/te_block_1 te_block_1.cpp

te_block_fixed: te_block_fixed.cpp
	mkdir -p $(BIN_DIR)
	g++ -O2 -Wall -lboost_program_options -o $(BIN_DIR)/te_block_fixed te_block_fixed.cpp

te_block: te_block.cpp
	mkdir -p $(BIN_DIR)
	g++ -O2 -Wall -lboost_program_options -o $(BIN_DIR)/te_block te_block.cpp

example: example.cpp
	mkdir -p $(BIN_DIR)
	g++ -O2 -Wall -o $(BIN_DIR)/example example.cpp

