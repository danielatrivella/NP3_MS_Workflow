SHELL = /bin/sh
CC = g++
CFLAGS = -O3 -std=c++11

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

.SUFFIXES: .cpp 

all:	$(OBJS)

.cpp.o:	$(OBJS) 
	$(CC) -c $(CFLAGS) -o $@ $<

clean-objs:
	@rm -f $(OBJS) 

clean-exe:
	@rm -f $(EXE) 

clean: clean-objs clean-exe
