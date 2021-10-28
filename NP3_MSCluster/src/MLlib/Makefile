SHELL = /bin/sh
CC = g++
CFLAGS = -O3 -std=c++11
LDFLAGS = -lm

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
COMMON_OBJS := $(patsubst %.cpp,%.o,$(wildcard ../Common/*.cpp))

EXE = mllib

.SUFFIXES: .cpp 

all: $(EXE)

$(EXE): $(OBJS) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(EXE) $(OBJS) $(COMMON_OBJS)

.cpp.o: 
	$(CC) -c $(CFLAGS) -o $@ $<

clean-objs:
	@rm -f $(OBJS) 

clean-exe:
	@rm -f $(EXE) 

clean: clean-objs clean-exe
