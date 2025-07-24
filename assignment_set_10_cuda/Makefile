# Compiler
NVCC      = nvcc
CXX       = g++

# Flags
NVCCFLAGS = -O3 -std=c++17 -Wno-deprecated-gpu-targets
CXXFLAGS  = -O3 -std=c++17

# Directories
SRC_DIR   = src
LIB_DIR   = lib
BUILD_DIR = build

# Files
TARGET    = main
CU_SOURCES = $(wildcard $(SRC_DIR)/*.cu) $(wildcard $(LIB_DIR)/*.cu)
CPP_SOURCES = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(LIB_DIR)/*.cpp)
HEADERS   = $(wildcard $(SRC_DIR)/*.hpp) $(wildcard $(LIB_DIR)/*.hpp) $(wildcard $(LIB_DIR)/*.cuh)

# Object files
SRC_CU_OBJS = $(patsubst $(SRC_DIR)/%.cu,$(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cu))
LIB_CU_OBJS = $(patsubst $(LIB_DIR)/%.cu,$(BUILD_DIR)/lib_%.o,$(wildcard $(LIB_DIR)/*.cu))
SRC_CPP_OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp))
LIB_CPP_OBJS = $(patsubst $(LIB_DIR)/%.cpp,$(BUILD_DIR)/lib_%.o,$(wildcard $(LIB_DIR)/*.cpp))

ALL_OBJS = $(SRC_CU_OBJS) $(LIB_CU_OBJS) $(SRC_CPP_OBJS) $(LIB_CPP_OBJS)

# Rules
all: $(BUILD_DIR) $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)
$(TARGET): $(ALL_OBJS)
	$(NVCC) $(NVCCFLAGS) -o $@ $^
# CUDA files from src/
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu $(HEADERS)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
# CUDA files from lib/
$(BUILD_DIR)/lib_%.o: $(LIB_DIR)/%.cu $(HEADERS)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# C++ files from src/
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# C++ files from lib/
$(BUILD_DIR)/lib_%.o: $(LIB_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean