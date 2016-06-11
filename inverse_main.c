#include <assert.h>
#include <stdio.h>
#include <stdlib.h>		// exit, EXIT_FAILURE

#include "inverse.h"

bool aliquot_inverse_cb(prime_pow_struct calc_list[], pow_idx_type pow_count) {
	//for (pow_idx_type i=0; i<pow_count; ++i)
	//	printf(PRI_NUM_TYPE "^" PRI_EXP_TYPE " ", calc_list[i].prime, calc_list[i].exp);
	//printf("\n");
	
	num_type value = 1;
	for (pow_idx_type i=0; i<pow_count; ++i) {
		// XXX req_aliquot_sum < 4294967312 = 2^32+16
		// for sum = 4294967312: n = 4294967311^2 > 2^64
		num_type pow = calc_pow(calc_list[i].prime, calc_list[i].exp);
		assert(value <= NUM_TYPE_MAX / pow);
		value *= pow;
	}
	printf(PRI_NUM_TYPE "\n", value);
	
	return false;
}

void run(const char *primes_filename) {
	aliquot_inverse_init(primes_filename);
	aliquot_inverse_set_callback(aliquot_inverse_cb);
	
	num_type req_aliquot_sum;
	int scanf_res = scanf(SCN_NUM_TYPE, &req_aliquot_sum);
	if (scanf_res != 1) {
		fprintf(stderr, "Wrong input\n");
		exit(EXIT_FAILURE);
	}
	aliquot_inverse_run(req_aliquot_sum);
	
	aliquot_inverse_terminate();
}

int main(int argc, char *argv[]) {
	if (argc != 2) {
		fprintf(stderr, "Too few arguments\n");
		exit(EXIT_FAILURE);
	}
	run(argv[1]);
	return 0;
}

