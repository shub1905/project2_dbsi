all:
	gcc -std=c99 -pedantic -msse4.2 -msse4a -o build main.c p2random.c tree.c
