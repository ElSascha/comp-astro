# Compiler
NVCC      = nvcc
CXX       = g++

# Flags
NVCCFLAGS = -O3 -std=c++17
CXXFLAGS  = -O3 -std=c++17

# Directories
SRC_DIR   = src
LIB_DIR   = lib
BUILD_DIR = build

# Dateien
TARGET    = main
CU_SOURCES = $(wildcard $(SRC_DIR)/*.cu)
CPP_SOURCES = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(LIB_DIR)/*.cpp)
HEADERS   = $(wildcard $(SRC_DIR)/*.hpp) $(wildcard $(LIB_DIR)/*.hpp)

# Object-Dateien
CU_OBJS   = $(CU_SOURCES:$(SRC_DIR)/%.cu=$(BUILD_DIR)/%.o)
CPP_OBJS  = $(CPP_SOURCES:%.cpp=$(BUILD_DIR)/%.o)

# Regeln
all: $(BUILD_DIR) $(TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TARGET): $(CU_OBJS) $(CPP_OBJS)
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu $(HEADERS)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(LIB_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)