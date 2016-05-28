#include <assert.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <fenv.h>

#include "typedefs.h"
#include "primes.h"

struct primes_array_struct_ {
	num_type *primes;
	size_t count;
};

primes_array_struct default_primes_array;

#define primes_buffer_size (1024*1024)
num_type primes_buffer[primes_buffer_size];

static size_t fill_primes(num_type primes[], size_t primes_size, num_type max_num) {
	if (primes_size == 0 || max_num < 2) return 0;
	int default_rounding_direction = fegetround();
	fesetround(0);
	primes[0] = 2;
	size_t primes_count = 1;
	num_type n = 3;
	num_type n_sqrt = 1;
	while (primes_count < primes_size && n <= max_num) {
		bool is_prime = true;
		for (size_t i=0; primes[i]<=n_sqrt; ++i) {
			if (!(n % primes[i])) {
				is_prime = false;
				break;
			}
		}
		assert(primes_count < primes_size);
		if (is_prime) primes[primes_count++] = n;
		n += 2;
		n_sqrt = round_sqrt(n);
	}
	fesetround(default_rounding_direction);
	return primes_count;
}

primes_array_struct *primes_construct() {
	size_t in_count = 1024*32; // TODO
	
	assert(in_count >= MAX_POW_COUNT);
	assert(in_count <= primes_buffer_size);
	size_t primes_count = fill_primes(primes_buffer, in_count, UINT64_MAX);
	assert(primes_count == in_count);
	
	default_primes_array.primes = primes_buffer;
	default_primes_array.count = primes_count;
	return &default_primes_array;
}

void primes_destruct(primes_array_struct *primes_array) {
	(void)primes_array;
}

num_type *primes_get_array(primes_array_struct *primes_array) {
	return primes_array->primes;
}

size_t primes_get_count(primes_array_struct *primes_array) {
	return primes_array->count;
}

