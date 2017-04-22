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

void primearr_construct(primearr_type *obj, const char *fpath, const prime_type *min_array, size_t min_size) {
	assert(obj != NULL);
	assert(min_array != NULL || min_size == 0);
	if (fpath != NULL) {
		int fd = open(fpath, O_RDONLY);
		if (fd < 0) {perror("open"); exit(EXIT_FAILURE);}
		struct stat file_stat;
		if (fstat(fd, &file_stat)) {perror("fstat"); exit(EXIT_FAILURE);}
		off_t file_size = file_stat.st_size;
		if (file_size % sizeof(prime_type) != 0) {fprintf(stderr, "Wrong file size\n"); exit(EXIT_FAILURE);}
		size_t file_arr_size = file_size / sizeof(prime_type);
		if (file_arr_size > min_size) {
			void *data = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
			obj->array = (const prime_type *)data;
			obj->size = file_arr_size;
			obj->is_mmapped = true;
		} else {
			obj->array = min_array;
			obj->size = min_size;
			obj->is_mmapped = false;
		}
		if (close(fd)) {perror("close"); exit(EXIT_FAILURE);}
	} else {
		obj->array = min_array;
		obj->size = min_size;
		obj->is_mmapped = false;
	}
}

void primearr_destruct(primearr_type *obj) {
	assert(obj != NULL);
	if (! obj->is_mmapped) return;
	if (munmap((void *)(obj->array), obj->size*sizeof(prime_type))) {perror("munmap"); exit(EXIT_FAILURE);}
}

const prime_type *primearr_get_array(const primearr_type *obj) {
	assert(obj != NULL);
	return obj->array;
}

size_t primearr_get_size(const primearr_type *obj) {
	assert(obj != NULL);
	return obj->size;
}
