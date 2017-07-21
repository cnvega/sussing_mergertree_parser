CC = h5cc
CCFLAGS = -std=c99 -Wno-unused-result -lm
CCFLAGS += -O2 

all: do_stats parse_mergertree

.PHONY: clean

%.o: %.c
	$(CC) -c $< -o $@ $(CCFLAGS)

do_stats: do_stats.o mergerTree.o
	$(CC) do_stats.o mergerTree.o -o do_stats $(CCFLAGS)

parse_mergertree: parse_mergertree.o mergerTree.o
	$(CC) parse_mergertree.o mergerTree.o -o parse_mergertree $(CCFLAGS)

clean:
	rm -f do_stats parse_mergertree *.o

