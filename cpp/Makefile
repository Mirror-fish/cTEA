# Makefile for cTEA C++ components
# Optimized for speed: single-pass scanning, in-memory data structures
# Avoids xTEA's chromosome-splitting bottleneck

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -fPIC -pthread

# For Apple Silicon macOS (Homebrew)
HTSLIB_DIR = /opt/homebrew/lib
HTSLIB_INC = /opt/homebrew/include

# For Intel macOS with Homebrew htslib:
# HTSLIB_DIR = /usr/local/lib
# HTSLIB_INC = /usr/local/include

# For Linux:
# HTSLIB_DIR = /usr/lib/x86_64-linux-gnu
# HTSLIB_INC = /usr/include/htslib

INCLUDES = -I$(HTSLIB_INC) -I.
LIBS = -L$(HTSLIB_DIR) -lhts -lz -lbz2 -llzma -lpthread

# For building shared/static library, define BUILD_LIB to exclude main functions
CXXFLAGS_LIB = $(CXXFLAGS) -DBUILD_LIB

# pybind11 configuration (use virtual environment's python)
PYBIND11_INC = $(shell /Users/jzeng/Desktop/cTEA/venv/bin/python3 -c "import pybind11; print(pybind11.get_include())" 2>/dev/null || echo "")
PYTHON_INC = /opt/homebrew/opt/python@3.14/Frameworks/Python.framework/Versions/3.14/include/python3.14
PYBIND11_CFLAGS = $(if $(PYBIND11_INC),-I$(PYBIND11_INC) -I$(PYTHON_INC),)

# Source files
SRCS = bam_processor.cpp extract_reads.cpp kmer_filter.cpp candidate_detector.cpp
OBJS = $(SRCS:.cpp=.o)
BINDINGS_OBJ = bindings.o

# Targets
ALL = libcTEA.a libcTEA.so libcTEA_pybind11.so extract_reads kmer_filter candidate_detector

all: $(ALL)

# Static library
libcTEA.a: $(OBJS)
	ar rcs $@ $^

# Shared library (for Python ctypes)
libcTEA.so: $(OBJS)
	$(CXX) -shared -o $@ $^ $(LIBS)

# pybind11 module (for direct Python import)
PYTHON_FRAMEWORK_PATH = /opt/homebrew/opt/python@3.14/Frameworks

libcTEA_pybind11.so: $(BINDINGS_OBJ) $(OBJS)
	$(CXX) -shared -o $@ $^ $(LIBS) -L$(HTSLIB_DIR) -lhts -lz -lbz2 -llzma -lpthread -F$(PYTHON_FRAMEWORK_PATH) -framework Python

# Object file for pybind11 bindings
bindings.o: bindings.cpp bam_processor.h candidate_detector.h
	$(CXX) $(CXXFLAGS) $(PYBIND11_CFLAGS) $(INCLUDES) -c $< -o $@

# Individual executables (for testing) - without BUILD_LIB
extract_reads: extract_reads_main.o bam_processor_main.o
	$(CXX) -o $@ $^ $(LIBS)

kmer_filter: kmer_filter_main.o
	$(CXX) -o $@ $^ $(LIBS)

candidate_detector: candidate_detector_main.o bam_processor_main.o
	$(CXX) -o $@ $^ $(LIBS)

# Object files for library (without main) - use BUILD_LIB
bam_processor.o: bam_processor.cpp bam_processor.h
	$(CXX) $(CXXFLAGS_LIB) $(INCLUDES) -c $< -o $@

extract_reads.o: extract_reads.cpp bam_processor.h
	$(CXX) $(CXXFLAGS_LIB) $(INCLUDES) -c $< -o $@

kmer_filter.o: kmer_filter.cpp
	$(CXX) $(CXXFLAGS_LIB) $(INCLUDES) -c $< -o $@

candidate_detector.o: candidate_detector.cpp bam_processor.h
	$(CXX) $(CXXFLAGS_LIB) $(INCLUDES) -c $< -o $@

# Object files for executables (with main) - no BUILD_LIB
bam_processor_main.o: bam_processor.cpp bam_processor.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

extract_reads_main.o: extract_reads.cpp bam_processor.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

kmer_filter_main.o: kmer_filter.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

candidate_detector_main.o: candidate_detector.cpp bam_processor.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(BINDINGS_OBJ) *_main.o $(ALL) libcTEA_pybind11.so

# Install (optional)
install: libcTEA.so
	cp libcTEA.so $(HTSLIB_DIR)/
	cp *.h /usr/local/include/cTEA/

.PHONY: all clean install