# aliquot
Calculate direct and inverse aliquot sum function.
[Aliquot sum](https://en.wikipedia.org/wiki/Aliquot_sum)
of a positive integer n is sum of all divisors of n other than n itself.

## inverse program
For given number `inverse` program finds all numbers which aliquot sum is equal to given number.
Result numbers are not in ascending order.

### Compilation
`make` or `mkdir -p build; make inverse`

### Usage
```
inverse [primes_file]
```
`primes_file` -- precomputed file of some first primes.
It has to contain some first prime numbers in 32-bit representation in ascending order starting from 2.
Amount of first 32-bit primes in file can be between 0 and 203280221
(203280221-st prime number 4294967291 = 0xFFFFFFFB is the largest 32-bit prime number).

Program reads input number from `STDIN` and prints result numbers to `STDOUT`.
Input number must be a positive integer greater than 1.

Example:
```
$ echo '9' | inverse
15
$ echo '8' | inverse primes.bin
49
10
```

## aliquot program
`aliquot` calculates aliquot sum of given number.

### Compilation
`mkdir -p build; make aliquot`

### Usage
Input and output formats and run arguments are similar to previous one.

Example:
```
$ echo '49' | aliquot
8
$ echo '10' | aliquot primes.bin
8
```
