#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "p2random.h"
#include "tree.h"

#define BILLION 1E9

int main(int argc, char* argv[]) {
        // parsing arguments
        assert(argc > 3);
        clock_t gen_start, gen_fin, simd_start, simd_fin, hard_start, hard_fin;
        double accum;
        size_t num_keys = strtoull(argv[1], NULL, 0);
        size_t num_probes = strtoull(argv[2], NULL, 0);
        size_t num_levels = (size_t) argc - 3;
        size_t* fanout = malloc(sizeof(size_t) * num_levels);
        assert(fanout != NULL);
        for (size_t i = 0; i < num_levels; ++i) {
                fanout[i] = strtoull(argv[i + 3], NULL, 0);
                assert(fanout[i] >= 2 && fanout[i] <= 17);
        }

        // building the tree index
        rand32_t* gen = rand32_init((uint32_t) time(NULL));
        assert(gen != NULL);
        int32_t* delimiter = generate_sorted_unique(num_keys, gen);
        assert(delimiter != NULL);
        Tree* tree = build_index(num_levels, fanout, num_keys, delimiter);
        free(delimiter);
        free(fanout);
        if (tree == NULL) {
                free(gen);
                exit(EXIT_FAILURE);
        }

        // generate probes
        int32_t* probe = generate(num_probes, gen);
        assert(probe != NULL);
        free(gen);
        uint32_t* result = malloc(sizeof(uint32_t) * num_probes);
        uint32_t* result_simd = malloc(sizeof(uint32_t) * num_probes);
        assert(result != NULL);
        assert(result_simd != NULL);

        // print_tree(tree);
        // printf("Probes : ");
        // for (int p=0; p<num_probes; p++)
        //         printf("%d ", probe[p]);
        // printf("\n");

        // perform index probing (Phase 2)
        // clock_gettime(CLOCK_REALTIME, &gen_start);
        gen_start = clock();
        for (size_t i = 0; i < num_probes; ++i) {
                result[i] = probe_index(tree, probe[i]);
        }
        gen_fin = clock();
        // clock_gettime(CLOCK_REALTIME, &gen_fin);

        // clock_gettime(CLOCK_REALTIME, &simd_start);
        simd_start = clock();
        for (size_t i = 0; i < num_probes; ++i) {
                result_simd[i] = probe_index_simd(tree, probe[i]);
        }
        simd_fin = clock();
        // clock_gettime(CLOCK_REALTIME, &simd_fin);

        // clock_gettime(CLOCK_REALTIME, &hard_start);
        if (tree->node_capacity[0] == 8 && tree->node_capacity[1] == 4 && tree->node_capacity[2] == 8 && num_probes%4==0) {
                hard_start = clock();
                uint32_t * result_simd_hardcode = probe_index_hardcode(tree, probe, num_probes);
                hard_fin = clock();
                accum = (double)(hard_fin - hard_start) / CLOCKS_PER_SEC;
                printf( "Hardcode Execution Time: %lf sec\n", accum );                
        }
        // clock_gettime(CLOCK_REALTIME, &hard_fin);

        // output results
        // for (size_t i = 0; i < num_probes; ++i) {
        //         fprintf(stdout, "%d %u %u %u\n", probe[i], result[i], result_simd[i], result_simd_hardcode[i]);
        // }

        // Calculate time it took
        // double accum = ( gen_fin.tv_sec - gen_start.tv_sec ) + ( gen_fin.tv_nsec - gen_start.tv_nsec ) / BILLION;
        accum = (double)(gen_fin - gen_start) / CLOCKS_PER_SEC;
        printf( "Regular Execution Time: %lf sec\n", accum );

        // accum = ( simd_fin.tv_sec - simd_start.tv_sec ) + ( simd_fin.tv_nsec - simd_start.tv_nsec ) / BILLION;
        accum = (double)(simd_fin - simd_start) / CLOCKS_PER_SEC;
        printf( "SIMD Execution Time: %lf sec\n", accum );

        // accum = ( hard_fin.tv_sec - hard_start.tv_sec ) + ( hard_fin.tv_nsec - hard_start.tv_nsec ) / BILLION;
        // cleanup and exit
        free(result);
        free(probe);
        cleanup_index(tree);
        return EXIT_SUCCESS;
}