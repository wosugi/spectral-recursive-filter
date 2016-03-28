# This is a -*- Makefile -*- 
# Warning: this file contains tabs which cannot be converted to spaces

CXX	=	$(CROSS)g++
LD	=	$(CXX)

# CFLAGS_OPT is normally in the config$(CONFIG).mk file,
# but here's a reasonable default value
# add -ftree-vectorizer-verbose=3 to get information about vectorized loops
CFLAGS_OPT= -O3 -Wall
#CFLAGS_OPENMP=-fopenmp

# the flag for .cpp files
#CXXFLAGS = -std=c++11

OPENCV_CPPFLAGS = $(shell pkg-config opencv --cflags)
OPENCV_LDFLAGS = $(shell pkg-config opencv --libs-only-L)
OPENCV_LIBS = $(shell pkg-config opencv --libs-only-l) $(shell pkg-config opencv --libs-only-other)

LIBTOOL = libtool

-include config$(CONFIG).mk

CPPFLAGS= $(OPENCV_CPPFLAGS)
LIBS = $(OPENCV_LIBS)


CFLAGS_EXTRA=$(CFLAGS_OPENMP)
CXXFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA)
CFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA)
LDFLAGS=$(CFLAGS_OPT) $(CFLAGS_EXTRA) $(LDFLAGS_EXTRA) $(OPENCV_LDFLAGS)

# for gcov profiling add:
#-fprofile-arcs -ftest-coverage

COMPILE.cpp=$(CXX) -c $(CXXFLAGS) $(CPPFLAGS)

PROGRAMS = srf

.PHONY: all

all: $(PROGRAMS)

# for each program, give list of sources and required libs
srf_SOURCES_CPP = \
	main.cpp
srf_HEADERS = \
	convolution.hpp \
	spectral_recursive_filter.hpp
srf_OBJS = $(srf_SOURCES_CPP:.cpp=.o)
srf_LIBS = $(LIBS)


# all sources and headers
SRCS_CPP = \
	$(srf_SOURCES_CPP)

HEADERS = \
	$(srf_HEADERS)

# link command for each program
srf: $(srf_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(srf_LIBS) $(LDADD)

.SUFFIXES: .o .cpp

## gcc-only version:

%.o : %.cpp
	$(COMPILE.cpp) -MD -o $@ $<
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

## general version:

# %.o : %.cpp
# 	@$(MAKEDEPEND); \
# 	    cp $*.d $*.P; \
# 	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
# 		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
# 	    rm -f $*.d
# 	$(COMPILE.cpp) -o $@ $<

.PHONY: clean distclean
clean:
	-rm -f $(PROGRAMS) *.o  *~
	-rm -rf *.dSYM

distclean: clean
	-rm -f $(SRCS_CPP:.cpp=.P)

count:
	 wc -l $(SRCS_CPP) $(HEADERS)

-include $(SRCS_CPP:.cpp=.P)
