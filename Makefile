all:
	gcc -std=gnu99 -pedantic -msse4.2 -msse4a -O2 -o build main.c p2random.c tree.c -lrt

# .PHONY: all
# all: clean build

# .PHONY: clean
# clean:
# 	rm -f *.o build
