#ifndef ALIQUOT_H
#define ALIQUOT_H

#include "typedefs.h"

void aliquot_init(const char *primes_filename);
void aliquot_terminate();
num_type calc_sigma(num_type num);
num_type calc_aliquot_sum(num_type num);

#endif/*ALIQUOT_H*/
