#ifndef PRIMES_H
#define PRIMES_H

#include <stddef.h>
#include "typedefs.h"

struct primes_array_struct_;
typedef struct primes_array_struct_ primes_array_struct;

primes_array_struct *primes_construct();
void primes_destruct(primes_array_struct *primes_array);
num_type *primes_get_array(primes_array_struct *primes_array);
size_t primes_get_count(primes_array_struct *primes_array);

#endif/*PRIMES_H*/

