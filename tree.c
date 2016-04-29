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
        size_t k;

        if (high == 4) {
            k = probe_index_simd_5(tree, probe_key, level, offset);
            result = result * (tree->node_capacity[level] + 1) + k;
        }
        else if (high == 8) {
            k = probe_index_simd_9(tree, probe_key, level, offset);
            result = result * (tree->node_capacity[level] + 1) + k;
        }
        else if (high == 16) {
            k = probe_index_simd_17(tree, probe_key, level, offset);
            result = result * (tree->node_capacity[level] + 1) + k;
        }
        else {
            printf("Error\n");
            return -1;
        }
    }
    return result;
}

inline uint32_t probe_index_simd_17(Tree* tree, int32_t probe_key, int32_t level, int32_t offset) {
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

inline uint32_t probe_index_simd_5(Tree* tree, int32_t probe_key, int32_t level, int32_t offset) {
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

inline uint32_t probe_index_simd_9(Tree* tree, int32_t probe_key, int32_t level, int32_t offset) {
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

uint32_t * probe_index_hardcode(Tree* tree, int32_t *probes, size_t num_probes){
    register __m128i root1 = _mm_load_si128((__m128i*) &tree->key_array[0][0]);
    register __m128i root2 = _mm_load_si128((__m128i*) &tree->key_array[0][4]);

    int mask1, mask2, mask3, mask4;
    int res1, res2, res3, res4;
    int res1_0 = 0, res2_0 = 0, res3_0 = 0, res4_0 = 0;

    register __m128i k1, k2, k3, k4;
    __m128i k;
    __m128i cmp_1, cmp_2, cmp_3, cmp_4;
    __m128i cmp_1_A, cmp_1_B, cmp_2_A, cmp_2_B;
    __m128i cmp_3_A, cmp_3_B, cmp_4_A, cmp_4_B;
    __m128i lvl1, lvl2, lvl3, lvl4;

    __m128i lvl_1_A, lvl_1_B, lvl_2_A, lvl_2_B;
    __m128i lvl_3_A, lvl_3_B, lvl_4_A, lvl_4_B;

    uint32_t* result = malloc(sizeof(uint32_t) * num_probes);

    for (int i=0; i<num_probes; i+=4) {
        // loading Keys
        k = _mm_load_si128((__m128i*) &probes[i]);
        k1 = _mm_shuffle_epi32(k, _MM_SHUFFLE(0,0,0,0));
        k2 = _mm_shuffle_epi32(k, _MM_SHUFFLE(1,1,1,1));
        k3 = _mm_shuffle_epi32(k, _MM_SHUFFLE(2,2,2,2));
        k4 = _mm_shuffle_epi32(k, _MM_SHUFFLE(3,3,3,3));

        // Scanning Level 0 capacity = 8
        cmp_1_A = _mm_cmpgt_epi32(k1, root1);
        cmp_1_B = _mm_cmpgt_epi32(k1, root2);

        cmp_2_A = _mm_cmpgt_epi32(k2, root1);
        cmp_2_B = _mm_cmpgt_epi32(k2, root2);

        cmp_3_A = _mm_cmpgt_epi32(k3, root1);
        cmp_3_B = _mm_cmpgt_epi32(k3, root2);                
        
        cmp_4_A = _mm_cmpgt_epi32(k4, root1);
        cmp_4_B = _mm_cmpgt_epi32(k4, root2);

        cmp_1 = _mm_packs_epi32(cmp_1_A, cmp_1_B);
        cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
        cmp_3 = _mm_packs_epi32(cmp_3_A, cmp_3_B);
        cmp_4 = _mm_packs_epi32(cmp_4_A, cmp_4_B);

        cmp_1 = _mm_packs_epi16(cmp_1, _mm_setzero_si128());
        cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
        cmp_3 = _mm_packs_epi16(cmp_3, _mm_setzero_si128());
        cmp_4 = _mm_packs_epi16(cmp_4, _mm_setzero_si128());

        mask1 = _mm_movemask_epi8(cmp_1);
        mask2 = _mm_movemask_epi8(cmp_2);
        mask3 = _mm_movemask_epi8(cmp_3);
        mask4 = _mm_movemask_epi8(cmp_4);

        res1_0 = _bit_scan_forward(mask1 ^ 0x1FFFF);
        res2_0 = _bit_scan_forward(mask2 ^ 0x1FFFF);
        res3_0 = _bit_scan_forward(mask3 ^ 0x1FFFF);
        res4_0 = _bit_scan_forward(mask4 ^ 0x1FFFF);

        res1 = res1_0;
        res2 = res2_0;
        res3 = res3_0;
        res4 = res4_0;

        // Scanning Level 1 capacity = 4
        lvl1 = _mm_load_si128((__m128i*) &tree->key_array[1][res1 << 2]);
        lvl2 = _mm_load_si128((__m128i*) &tree->key_array[1][res2 << 2]);
        lvl3 = _mm_load_si128((__m128i*) &tree->key_array[1][res3 << 2]);
        lvl4 = _mm_load_si128((__m128i*) &tree->key_array[1][res4 << 2]);

        cmp_1 = _mm_cmpgt_epi32(k1, lvl1);
        cmp_2 = _mm_cmpgt_epi32(k2, lvl2);
        cmp_3 = _mm_cmpgt_epi32(k3, lvl3);
        cmp_4 = _mm_cmpgt_epi32(k4, lvl4);

        mask1 = _mm_movemask_ps(_mm_castsi128_ps(cmp_1));
        mask2 = _mm_movemask_ps(_mm_castsi128_ps(cmp_2));
        mask3 = _mm_movemask_ps(_mm_castsi128_ps(cmp_3));
        mask4 = _mm_movemask_ps(_mm_castsi128_ps(cmp_4));

        res1_0 = _bit_scan_forward(mask1 ^ 0x1FF);
        res2_0 = _bit_scan_forward(mask2 ^ 0x1FF);
        res3_0 = _bit_scan_forward(mask3 ^ 0x1FF);
        res4_0 = _bit_scan_forward(mask4 ^ 0x1FF);

        res1 = res1 * 5 + res1_0;
        res2 = res2 * 5 + res2_0;
        res3 = res3 * 5 + res3_0;
        res4 = res4 * 5 + res4_0;

        //Scanning Level 2 capacity = 8
        lvl_1_A = _mm_load_si128((__m128i*) &tree->key_array[2][res1 << 3]);
        lvl_1_B = _mm_load_si128((__m128i*) &tree->key_array[2][(res1 << 3) + 4]);

        lvl_2_A = _mm_load_si128((__m128i*) &tree->key_array[2][res2 << 3]);
        lvl_2_B = _mm_load_si128((__m128i*) &tree->key_array[2][(res2 << 3) + 4]);

        lvl_3_A = _mm_load_si128((__m128i*) &tree->key_array[2][res3 << 3]);
        lvl_3_B = _mm_load_si128((__m128i*) &tree->key_array[2][(res3 << 3) + 4]);

        lvl_4_A = _mm_load_si128((__m128i*) &tree->key_array[2][res4 << 3]);
        lvl_4_B = _mm_load_si128((__m128i*) &tree->key_array[2][(res4 << 3) + 4]);

        cmp_1_A = _mm_cmpgt_epi32(k1, lvl_1_A);
        cmp_1_B = _mm_cmpgt_epi32(k1, lvl_1_B);

        cmp_2_A = _mm_cmpgt_epi32(k2, lvl_2_A);
        cmp_2_B = _mm_cmpgt_epi32(k2, lvl_2_B);

        cmp_3_A = _mm_cmpgt_epi32(k3, lvl_3_A);
        cmp_3_B = _mm_cmpgt_epi32(k3, lvl_3_B);

        cmp_4_A = _mm_cmpgt_epi32(k4, lvl_4_A);
        cmp_4_B = _mm_cmpgt_epi32(k4, lvl_4_B);

        cmp_1 = _mm_packs_epi32(cmp_1_A, cmp_1_B);
        cmp_2 = _mm_packs_epi32(cmp_2_A, cmp_2_B);
        cmp_3 = _mm_packs_epi32(cmp_3_A, cmp_3_B);
        cmp_4 = _mm_packs_epi32(cmp_4_A, cmp_4_B);

        cmp_1 = _mm_packs_epi16(cmp_1, _mm_setzero_si128());
        cmp_2 = _mm_packs_epi16(cmp_2, _mm_setzero_si128());
        cmp_3 = _mm_packs_epi16(cmp_3, _mm_setzero_si128());
        cmp_4 = _mm_packs_epi16(cmp_4, _mm_setzero_si128());

        mask1 = _mm_movemask_epi8(cmp_1);
        mask2 = _mm_movemask_epi8(cmp_2);
        mask3 = _mm_movemask_epi8(cmp_3);
        mask4 = _mm_movemask_epi8(cmp_4);

        res1_0 = _bit_scan_forward(mask1 ^ 0x1FFFF);
        res2_0 = _bit_scan_forward(mask2 ^ 0x1FFFF);
        res3_0 = _bit_scan_forward(mask3 ^ 0x1FFFF);
        res4_0 = _bit_scan_forward(mask4 ^ 0x1FFFF);

        res1 = res1 * 9 + res1_0;
        res2 = res2 * 9 + res2_0;
        res3 = res3 * 9 + res3_0;
        res4 = res4 * 9 + res4_0;

        result[i] = res1;
        result[i+1] = res2;
        result[i+2] = res3;
        result[i+3] = res4;
    }
    return result;
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

double time_difference(struct timespec start, struct timespec finish) {
	double elapsed;
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	return elapsed;
}
