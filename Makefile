CC=gcc
LD=gcc
STRIP=strip -s
#STRIP=true
WARNINGS=-Wall -Wextra -pedantic
DEBUG=
#DEBUG=-g -ggdb
#COPTIM=-march=native -O3
COPTIM=-march=native -O2
#COPTIM=-O0
#DEFINES=-DNDEBUG
DEFINES=
INCLUDES=
CSTD=-std=c11
CFLAGS=$(WARNINGS) $(DEBUG) $(COPTIM) $(DEFINES) $(INCLUDES) $(CSTD) -pipe
LDOPTIM=-Wl,-O1 -Wl,--as-needed
#LDOPTIM=
LIBFILES=-lm
LDFLAGS=$(WARNINGS) $(DEBUG) $(LDOPTIM) $(LIBFILES)
SRC_DIR=.
BUILD_DIR=build

.PHONY: all clean

all: $(BUILD_DIR) inverse

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

aliquot: $(BUILD_DIR)/aliquot.o
	$(LD) -o $@ $^ $(LDFLAGS)
	$(STRIP) $@

$(BUILD_DIR)/aliquot.o: $(SRC_DIR)/aliquot.c Makefile
	$(CC) -o $@ $< -c $(CFLAGS)

inverse: $(BUILD_DIR)/inverse_main.o $(BUILD_DIR)/inverse.o $(BUILD_DIR)/primes.o
	$(LD) -o $@ $^ $(LDFLAGS)
	$(STRIP) $@

$(BUILD_DIR)/inverse_main.o: $(SRC_DIR)/inverse_main.c $(SRC_DIR)/inverse.h $(SRC_DIR)/typedefs.h Makefile
	$(CC) -o $@ $< -c $(CFLAGS)

$(BUILD_DIR)/inverse.o: $(SRC_DIR)/inverse.c $(SRC_DIR)/inverse.h $(SRC_DIR)/primes.h $(SRC_DIR)/typedefs.h Makefile
	$(CC) -o $@ $< -c $(CFLAGS)

$(BUILD_DIR)/primes.o: $(SRC_DIR)/primes.c $(SRC_DIR)/primes.h Makefile
	$(CC) -o $@ $< -c $(CFLAGS)

clean:
	rm -r $(BUILD_DIR)

