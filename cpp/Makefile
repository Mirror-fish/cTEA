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

# Source files
SRCS = bam_processor.cpp extract_reads.cpp kmer_filter.cpp candidate_detector.cpp
OBJS = $(SRCS:.cpp=.o)

# Targets
ALL = libcTEA.a libcTEA.so extract_reads kmer_filter candidate_detector

all: $(ALL)

# Static library
libcTEA.a: $(OBJS)
	ar rcs $@ $^

# Shared library (for Python ctypes/pybind11)
libcTEA.so: $(OBJS)
	$(CXX) -shared -o $@ $^ $(LIBS)

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
	rm -f $(OBJS) *_main.o $(ALL)

# Install (optional)
install: libcTEA.so
	cp libcTEA.so $(HTSLIB_DIR)/
	cp *.h /usr/local/include/cTEA/

.PHONY: all clean install