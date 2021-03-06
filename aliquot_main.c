#include <assert.h>
#include <stdio.h>
#include <stdlib.h>		// exit, EXIT_FAILURE, NULL

#include "aliquot.h"

void run_aliquot_sequence() {
	num_type n;
	int scanf_res = scanf(SCN_NUM_TYPE, &n);
	if (scanf_res != 1) {
		fprintf(stderr, "Wrong input\n");
		exit(EXIT_FAILURE);
	}
	printf(PRI_NUM_TYPE "\n", n);
	
	while (n != 1) {
		n = calc_aliquot_sum(n);
		printf(PRI_NUM_TYPE "\n", n);
	}
}

void run_aliquot_sum_table() {
	num_type num_last;
	int scanf_res = scanf(SCN_NUM_TYPE, &num_last);
	if (scanf_res != 1) {
		fprintf(stderr, "Wrong input\n");
		exit(EXIT_FAILURE);
	}

	for (num_type n=2; n<=num_last; ++n) {
		//printf(PRI_NUM_TYPE "\t", n);
		printf(PRI_NUM_TYPE "\n", calc_aliquot_sum(n));
	}
}

void run_aliquot_sum() {
	num_type n;
	int scanf_res = scanf(SCN_NUM_TYPE, &n);
	if (scanf_res != 1) {
		fprintf(stderr, "Wrong input\n");
		exit(EXIT_FAILURE);
	}
	if (!(n > 1)) {
		fprintf(stderr, "Wrong input: " PRI_NUM_TYPE "\n", n);
		exit(EXIT_FAILURE);
	}
	printf(PRI_NUM_TYPE "\n", calc_aliquot_sum(n));
}

void run(const char *primes_filename) {
	aliquot_init(primes_filename);
	
	run_aliquot_sum();
	
	aliquot_terminate();
}

int main(int argc, char *argv[]) {
	if (argc > 2) {
		fprintf(stderr, "Too many arguments\n");
		exit(EXIT_FAILURE);
	}
	run(argc > 1 ? argv[1] : NULL);
	return 0;
}
