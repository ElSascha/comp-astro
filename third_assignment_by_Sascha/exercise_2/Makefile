# Compiler
CXX = g++
# Compiler flags

CXXFLAGS = -std=c++11 -Wall -O2 $(GSL_CXXFLAGS)
LDFLAGS = $(GSL_LDFLAGS) -lgsl -lgslcblas -lm
# Source files: main.cpp and all .cpp files in lib
SRC = main.cpp $(wildcard lib/*.cpp)
# Build directory
BUILD_DIR = build
# Object files
OBJ = $(SRC:%.cpp=$(BUILD_DIR)/%.o)
# Executable name
EXEC = main
# Default target
all: $(EXEC)
# Link object files to create the executable
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
# Compile source files into object files
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@
# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(EXEC)
.PHONY: all clean
