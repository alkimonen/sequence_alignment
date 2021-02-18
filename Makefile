all: sequenceAlignment
sequenceAlignment: sequenceAlignment.c
	gcc -g -o sequenceAlignment sequenceAlignment.c
clean:
	rm -fr sequenceAlignment sequenceAlignment.o *~
