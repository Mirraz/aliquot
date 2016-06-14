#include <assert.h>
#include <stdio.h>
#include <stdlib.h>		// exit, EXIT_FAILURE
#include <inttypes.h>

#include "aliquot.h"

void run_aliquot_sequence() {
	num_type n;
	int scanf_res = scanf("%" SCNuFAST64, &n);
	assert(scanf_res == 1);
	printf("%" PRIuFAST64 "\n", n);
	
	while (n != 1) {
		n = calc_aliquot_sum(n);
		printf("%" PRIuFAST64 "\n", n);
	}
}

void run_aliquot_sum_table() {
	num_type n;
	for (n=2; n<=16777216; ++n) {
		//printf("%" PRIuFAST64 "\t", n);
		printf("%" PRIuFAST64 "\n", calc_aliquot_sum(n));
	}
}

void run_aliquot_sum() {
	num_type n;
	int scanf_res = scanf("%" SCNuFAST64, &n);
	assert(scanf_res == 1);
	
	printf("%" PRIuFAST64 "\n", calc_aliquot_sum(n));
}

void run(const char *primes_filename) {
	aliquot_init(primes_filename);
	
	run_aliquot_sum();
	
	aliquot_terminate();
}

int main(int argc, char *argv[]) {
	if (argc != 2) {
		fprintf(stderr, "Too few arguments\n");
		exit(EXIT_FAILURE);
	}
	run(argv[1]);
	return 0;
}
