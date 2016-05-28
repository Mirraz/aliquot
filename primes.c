#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "primes.h"

struct primes_array_struct_ {
	prime_type *primes;
	size_t count;
};

primes_array_struct *primes_construct() {
	const char *primes_filename = "primes32.bin";
	
	primes_array_struct *primes_array = malloc(sizeof(primes_array_struct));
	if (primes_array == NULL) {perror("malloc"); exit(EXIT_FAILURE);}
	
	int fd = open(primes_filename, O_RDONLY);
	if (fd < 0) {perror("open"); exit(EXIT_FAILURE);}
	struct stat file_stat;
	if (fstat(fd, &file_stat)) {perror("fstat"); exit(EXIT_FAILURE);}
	off_t file_size = file_stat.st_size;
	if (file_size == 0 || file_size % sizeof(prime_type) != 0) {fprintf(stderr, "Wrong file size\n"); exit(EXIT_FAILURE);}
	void *data = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (close(fd)) {perror("close"); exit(EXIT_FAILURE);}
	
	primes_array->primes = (prime_type *)data;
	primes_array->count  = file_size / sizeof(prime_type);
	return primes_array;
}

void primes_destruct(primes_array_struct *primes_array) {
	if (munmap(primes_array->primes, primes_array->count*sizeof(prime_type))) {
		perror("munmap");
		exit(EXIT_FAILURE);
	}
	free(primes_array);
}

prime_type *primes_get_array(primes_array_struct *primes_array) {
	return primes_array->primes;
}

size_t primes_get_count(primes_array_struct *primes_array) {
	return primes_array->count;
}

