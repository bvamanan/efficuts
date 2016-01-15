CPP = g++
CFLGA = -lm -o3

compressedcuts: compressedcuts.c compressedcuts.h
	${CPP} ${CFLGA} -o compressedcuts compressedcuts.c

all: compressedcuts

clean: 
	rm -f compressedcuts
