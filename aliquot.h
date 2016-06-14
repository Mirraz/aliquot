#ifndef ALIQUOT_H
#define ALIQUOT_H

typedef uint_fast64_t num_type;

void aliquot_init(const char *primes_filename);
void aliquot_terminate();
num_type calc_sigma(num_type num);
num_type calc_aliquot_sum(num_type num);

#endif/*ALIQUOT_H*/

