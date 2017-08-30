################################################################################
## BEGIN: USER DEFINED            
################################################################################
CXX=g++
CC=gcc
CXX_FLAGS =-fPIC -std=c++14
RELEASE_FLAG=-O3 -Wall -Wno-unused-local-typedefs
DEBUG_FLAG=-g3 -Wall -Wno-unused-local-typedefs -DUSE_DEBUG
DEFINES=
LDFLAGS=
################################################################################
## END: USER DEFINED
################################################################################
SRC_DIR=src
BIN_DIR=bin
EXT_DIR=ext
UFLIO_SRC_DIR=$(SRC_DIR)/ufl_io
TESTS_DIR=$(SRC_DIR)/tests
BUILD_DIR=build
LIB_NAME=usymqr
LIB_DIR=$(BUILD_DIR)/lib
LIB=$(LIB_DIR)/lib$(LIB_NAME).so
OUT_DIR=$(BIN_DIR) $(BUILD_DIR) $(LIB_DIR)
CPP_FILES=$(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES=$(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(CPP_FILES))
OBJ_UFLIO=$(BUILD_DIR)/mmio.o
SRC_UFLIO=$(UFLIO_SRC_DIR)/mmio.cpp
INCLUDES =-I$(SRC_DIR)
INCLUDES+=-I$(EXT_DIR)/eigen # Eigen
################################################################################
## BUILD INSTRUCTIONS
################################################################################
.PHONY: all
all: directories release
.PHONY: debug
debug: CXX_FLAGS += $(DEBUG_FLAG) 
debug: directories tests
.PHONY: release
release: CXX_FLAGS += $(RELEASE_FLAG)
release: directories tests
.PHONY: directories
MKDIR_P = mkdir -p
directories: ${OUT_DIR}
${OUT_DIR}:
	${MKDIR_P} ${OUT_DIR}
######################q##########################################################
## LIBRARY & OBJ FILES
################################################################################
## TODO need to test if $(OBJ_FILES) is empty, if yes then don't build lib
$(LIB): $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) -shared -o $@ $(OBJ_FILES) $(LIBS) $(LIBS_DIRS)

$(OBJ_FILES): $(CPP_FILES)
	$(CXX) $(CXX_FLAGS) $(DEFINES) $(INCLUDES) -o $@ \
		-c $(patsubst $(BUILD_DIR)/%.o,$(SRC_DIR)/%.cpp,$@) $(LIBS) $(LIBS_DIRS)

$(OBJ_UFLIO): $(SRC_UFLIO)
	$(CXX) $(CXX_FLAGS) $(DEFINES) $(SRC_UFLIO) -o $@ -c
################################################################################
## INDIVIDUAL TESTS
################################################################################
tests: \
	$(OBJ_UFLIO) \
	usymqr_test \
	random_matrices_test \
	ufl_test
usymqr_test: $(LIB) $(BIN_DIR)
	$(CXX) $(CXX_FLAGS) $(DEFINES) $(INCLUDES) $(TESTS_DIR)/$@.cpp \
		-o $(BIN_DIR)/$@ $(LIBS) -l$(LIB_NAME) $(LIBS_DIRS) -L$(LIB_DIR)
random_matrices_test: $(LIB) $(BIN_DIR)
	$(CXX) $(CXX_FLAGS) $(DEFINES) $(INCLUDES) $(TESTS_DIR)/$@.cpp \
		-o $(BIN_DIR)/$@ $(LIBS) -l$(LIB_NAME) $(LIBS_DIRS) -L$(LIB_DIR)
ufl_test: $(LIB) $(BIN_DIR) $(OBJ_UFLIO)
	$(CXX) $(CXX_FLAGS) $(DEFINES) $(INCLUDES) $(TESTS_DIR)/$@.cpp \
		-o $(BIN_DIR)/$@ $(LIBS) -l$(LIB_NAME) $(LIBS_DIRS) -L$(LIB_DIR) \
		   $(OBJ_UFLIO)
	       
################################################################################
## CLEAN
################################################################################
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
################################################################################
## RUN
################################################################################
run: 
	bin/usymqr_test
