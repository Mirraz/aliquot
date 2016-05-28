#ifndef PRIMES_H
#define PRIMES_H

#include <stdint.h>
#include <stddef.h>

struct primes_array_struct_;
typedef struct primes_array_struct_ primes_array_struct;

typedef uint32_t prime_type;

primes_array_struct *primes_construct();
void primes_destruct(primes_array_struct *primes_array);
prime_type *primes_get_array(primes_array_struct *primes_array);
size_t primes_get_count(primes_array_struct *primes_array);

#endif/*PRIMES_H*/

