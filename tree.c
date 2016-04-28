#include "tree.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int posix_memalign(void** memptr, size_t alignment, size_t size);
size_t alignment = 16;

Tree* build_index(size_t num_levels, size_t fanout[], size_t num_keys, int32_t key[]) {
    // return null pointer for invalid tree configuration
    size_t min_num_keys = 1;
    for (size_t i = 0; i < num_levels - 1; ++i) {
        min_num_keys *= fanout[i];
    }
    size_t max_num_keys = min_num_keys * fanout[num_levels - 1] - 1;
    if (num_keys < min_num_keys || num_keys > max_num_keys) {
        fprintf(stderr, "Error: incorrect number of keys, min %zu, max %zu\n", min_num_keys, max_num_keys);
        return NULL;
    }

    // initialize the tree index
    Tree* tree = malloc(sizeof(Tree));
    assert(tree != NULL);
    tree->num_levels = num_levels;
    tree->node_capacity = malloc(sizeof(size_t) * num_levels);
    assert(tree->node_capacity != NULL);
    for (size_t i = 0; i < num_levels; ++i) {
        tree->node_capacity[i] = fanout[i] - 1;
    }
    tree->key_array = malloc(sizeof(int32_t*) * num_levels);
    assert(tree->key_array != NULL);

    tree->key_count = malloc(sizeof(size_t) * num_levels);
    assert(tree->key_count != NULL);
    size_t* key_count = malloc(sizeof(size_t) * num_levels);
    assert(key_count != NULL);

    size_t* array_capacity = malloc(sizeof(size_t) * num_levels);
    assert(array_capacity != NULL);
    for (size_t i = 0; i < num_levels; ++i) {
        size_t size = sizeof(int32_t) * tree->node_capacity[i];         // allocate one node per level
        int error = posix_memalign((void**) &(tree->key_array[i]), alignment, size);
        assert(error == 0);
        key_count[i] = 0;
        tree->key_count[i] = 0;
        array_capacity[i] = tree->node_capacity[i];     // array_capacity[i] is always a multiple of node_capacity[i]
    }

    // insert sorted keys into index
    for (size_t i = 1; i < num_keys; ++i) {
        assert(key[i - 1] < key[i]);
    }
    for (size_t i = 0; i < num_keys; ++i) {
        size_t level = num_levels - 1;
        while (key_count[level] == array_capacity[level])
            level -= 1;
        tree->key_array[level][key_count[level]] = key[i];
        key_count[level] += 1;
        tree->key_count[level] += 1;
        while (level < num_levels - 1) {
            level += 1;
            size_t new_capacity = array_capacity[level] + tree->node_capacity[level];
            size_t size = sizeof(int32_t) * new_capacity;           // allocate one more node
            int32_t* new_array = NULL;
            int error = posix_memalign((void**) &new_array, alignment, size);
            assert(error == 0);
            memcpy(new_array, tree->key_array[level], sizeof(int32_t) * key_count[level]);
            free(tree->key_array[level]);
            tree->key_array[level] = new_array;
            array_capacity[level] = new_capacity;
        }
    }

    // pad with INT32_MAXs
    for (size_t i = 0; i < num_levels; ++i) {
        for (size_t j = key_count[i]; j < array_capacity[i]; ++j)
            tree->key_array[i][j] = INT32_MAX;
        key_count[i] = array_capacity[i];
        tree->key_count[i] = array_capacity[i];
    }

    // print the tree
    // for (size_t i = 0; i < num_levels; ++i) {
    //         printf("Level %zu:", i);
    //         for (size_t j = 0; j < key_count[i]; ++j)
    //                 printf(" %d", tree->key_array[i][j]);
    //         printf("\n");
    // }

    free(array_capacity);
    free(key_count);
    return tree;
}

uint32_t probe_index(Tree* tree, int32_t probe_key) {
    size_t result = 0;
    for (size_t level = 0; level < tree->num_levels; ++level) {
        size_t offset = result * tree->node_capacity[level];
        size_t low = 0;
        size_t high = tree->node_capacity[level];
        while (low != high) {
            size_t mid = (low + high) / 2;
            if (tree->key_array[level][mid + offset] < probe_key)
                low = mid + 1;
            else
                high = mid;
        }
        size_t k = low;       // should go to child k
        result = result * (tree->node_capacity[level] + 1) + k;
    }
    return (uint32_t) result;
}

uint32_t probe_index_simd(Tree* tree, int32_t probe_key) {
    int result = 0;
    for (size_t level = 0; level < tree->num_levels; ++level) {
        size_t offset = result * tree->node_capacity[level];
        size_t low = 0;
        size_t high = tree->node_capacity[level];

        if (high == 4) {
            size_t k = probe_index_simd_5(tree, probe_key, level, offset);
            result = result * (tree->node_capacity[level] + 1) + k;
        }
        else if (high == 8) {
            size_t k = probe_index_simd_9(tree, probe_key, level, offset);
            result = result * (tree->node_capacity[level] + 1) + k;
        }
        else if (high == 16) {
            size_t k = probe_index_simd_17(tree, probe_key, level, offset);
            result = result * (tree->node_capacity[level] + 1) + k;
        }
        else {
            printf("Error\n");
            return -1;
        }
    }
}

uint32_t probe_index_simd_17(Tree* tree, int32_t probe_key, int32_t level, int32_t offset) {
    __m128i key = _mm_cvtsi32_si128(probe_key);
    key = _mm_shuffle_epi32(key, 0);

    __m128i del_ABCD = _mm_load_si128((__m128i*) &tree->key_array[level][offset + 0]);
    __m128i del_EFGH = _mm_load_si128((__m128i*) &tree->key_array[level][offset + 4]);
    __m128i del_IJKL = _mm_load_si128((__m128i*) &tree->key_array[level][offset + 8]);
    __m128i del_MNOP = _mm_load_si128((__m128i*) &tree->key_array[level][offset + 12]);


    // compare with 16 delimiters stored in 4 registers
    __m128i cmp_ABCD = _mm_cmpgt_epi32(key, del_ABCD);
    __m128i cmp_EFGH = _mm_cmpgt_epi32(key, del_EFGH);
    __m128i cmp_IJKL = _mm_cmpgt_epi32(key, del_IJKL);
    __m128i cmp_MNOP = _mm_cmpgt_epi32(key, del_MNOP);
    // pack results to 16-bytes in a single SIMD register
    __m128i cmp_A_to_H = _mm_packs_epi32(cmp_ABCD, cmp_EFGH);
    __m128i cmp_I_to_P = _mm_packs_epi32(cmp_IJKL, cmp_MNOP);
    __m128i cmp_A_to_P = _mm_packs_epi16(cmp_A_to_H, cmp_I_to_P);
    // extract the mask the least significant bit
    int mask = _mm_movemask_epi8(cmp_A_to_P);
    uint32_t res = _bit_scan_forward(mask ^ 0xFFFFFFFF);   
    return res; 
}

uint32_t probe_index_simd_5(Tree* tree, int32_t probe_key, int32_t level, int32_t offset) {
    // access level 1 (non-root) of the index (5-way)
    __m128i key = _mm_cvtsi32_si128(probe_key);
    key = _mm_shuffle_epi32(key, 0);

    __m128i lvl = _mm_load_si128((__m128i*) &tree->key_array[level][offset]);
    __m128i cmp = _mm_cmpgt_epi32(key, lvl);
    int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));

    // printf("Check : %d\n", tree->key_array[level][offset]);
    uint32_t res = _bit_scan_forward(mask ^ 0x1FF);
    // r_1 += (r_0 << 2) + r_0;
    return res;
}

uint32_t probe_index_simd_9(Tree* tree, int32_t probe_key, int32_t level, int32_t offset) {
    // access level 1 (non-root) of the index (5-way)
    __m128i key = _mm_cvtsi32_si128(probe_key);
    key = _mm_shuffle_epi32(key, 0);

    __m128i lvl_2_A = _mm_load_si128((__m128i*) &tree->key_array[level][offset]);
    __m128i lvl_2_B = _mm_load_si128((__m128i*) &tree->key_array[level][offset + 4]);

    __m128i cmp_2_A = _mm_cmpgt_epi32(key, lvl_2_A);
    __m128i cmp_2_B = _mm_cmpgt_epi32(key, lvl_2_B);
    
    __m128i cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
    cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
    int mask = _mm_movemask_epi8(cmp_2);
    int res = _bit_scan_forward(mask ^ 0x1FFFF);
    // r_2 += (r_1 << 3) + r_1;
    return res;
}

void print_tree(Tree* tree) {
    for (int i=0; i<tree->num_levels; i++) {
        // printf("Level %d %zu ", i, tree->key_count[i]);
        for (size_t j=0; j<tree->key_count[i]; j++) {
            printf("(%zu, %d) ", j, tree->key_array[i][j]);
        }
        printf("\n");
    }
}

void cleanup_index(Tree* tree) {
    free(tree->node_capacity);
    for (size_t i = 0; i < tree->num_levels; ++i)
        free(tree->key_array[i]);
    free(tree->key_array);
    free(tree);
}
