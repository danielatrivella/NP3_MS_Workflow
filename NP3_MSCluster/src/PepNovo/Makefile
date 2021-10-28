SHELL = /bin/sh
CC = g++
CFLAGS = -O3 -std=c++11
LDFLAGS = -lm

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
ALL_MLLIB_OBJS := $(patsubst %.cpp,%.o,$(wildcard ../MLlib/*.cpp))
MLLIB_OBJS := $(filter-out ../MLlib/mllib.o,$(ALL_MLLIB_OBJS))
COMMON_OBJS := $(patsubst %.cpp,%.o,$(wildcard ../Common/*.cpp))

.SUFFIXES: .cpp 

all: $(OBJS)

.cpp.o: 
	$(CC) -c $(CFLAGS) -o $@ $<

clean-objs:
	@rm -f $(OBJS) 

clean-exe:
	@rm -f $(EXE) 

clean: clean-objs clean-exe

print_vars:
	echo "OBJS: " $(OBJS)
	echo "ALL_MLLIB_OBJS: " $(ALL_MLLIB_OBJS)
	echo "MLLIB_OBJS: " $(MLLIB_OBJS)
	echo "COMMIN_OBJS: " $(COMMON_OBJS)

