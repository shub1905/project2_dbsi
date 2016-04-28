#ifndef TREE_H_
#define TREE_H_

#include <stdint.h>
#include <stdlib.h>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <nmmintrin.h>
#include <ammintrin.h>
#include <x86intrin.h>

typedef struct {
        size_t num_levels;
        size_t* node_capacity;
        int32_t** key_array;
} Tree;

Tree* build_index(size_t num_levels, size_t fanout[], size_t num_keys, int32_t key[]);
uint32_t probe_index(Tree* tree, int32_t probe_key);
uint32_t probe_index_simd(Tree* tree, int32_t probe_key, int32_t level, int32_t offset);
void cleanup_index(Tree* tree);

#endif
