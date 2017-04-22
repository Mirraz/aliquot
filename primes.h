#ifndef PRIMES_H
#define PRIMES_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

typedef uint32_t prime_type;
typedef struct {
	const prime_type *array;
	size_t size;
	bool is_mmapped;
} primearr_type;

void primearr_construct(primearr_type *obj, const char *fpath, const prime_type *min_array, size_t min_size);
void primearr_destruct(primearr_type *obj);
const prime_type *primearr_get_array(const primearr_type *obj);
size_t primearr_get_size(const primearr_type *obj);

#endif/*PRIMES_H*/

