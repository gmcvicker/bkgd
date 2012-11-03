#include <math.h>
#include <glib.h>
#include <limits.h>
#include <stdio.h>


/* numerical utilities */

/**
 * Returns the nbit most significant bits of the provided unsigned
 * integer.
 *
 * For example numer_sigbit(43, 2) will give the two most significant
 * bits of 000101011 or  32 = (00010000)
 *
 */
unsigned int numer_sigbit(const unsigned int num, const int nbit) {
  unsigned int a,max_bit;
  int i;

  if(nbit == 0) {
    return 0;
  }

  /* find number of bits in num */
  max_bit = UINT_MAX >> 1;

  if(nbit >= max_bit) {
    return num;
  }

  a = 1 << (nbit-1);

  i = nbit-1;
  /* find most significant bit in num */
  while(num >= a && a < max_bit) {
    a <<= 1;
    i++;
  }
  
  if(i-nbit > 0) {
    /* AND with a mask has ones in significant bits that we
     * want to keep, but 0s in lower order bits
     */
    return num & (UINT_MAX << (i-nbit));
  }

  /* there were less than or equal to nbit significant bits */
  return num;
}


/**
 * Same as numer_sigbit() but for unsigned long integers instead of
 * unsigned integers.
 *
 */
unsigned long numer_sigbitl(const unsigned long num, const int nbit) {
  unsigned long a,max_bit;
  int i;

  if(nbit == 0) {
    return 0;
  }

  /* find number of bits in num */
  max_bit = ULONG_MAX >> 1;

  /* find most significant bit in num */
  i = 0;
  a = 1;
  while(num >= a && a < max_bit) {
    a <<= 1;
    i++;
  }
  
  if(i-nbit > 0) {
    /* AND with a mask has ones in significant bits that we
     * want to keep, but 0s in lower order bits
     */
    return num & (ULONG_MAX << (i-nbit));
  }

  /* there were less than or equal to nbit significant bits */
  return num;
}



/*
 * Returns an approximation the the natural logarithm of the factorial
 * ln(n!) = ln(n) +  ln(n-1) + ln(n-2) + ... + ln(2) + ln(1)
 *
 * Uses Stirling's approximation for n >= 15.  
 * (see: http://mathworld.wolfram.com/StirlingsApproximation.html)
 *
 * ln(n!) ~= n * ln(n) - n
 * For speed and accuracy, returns pre-computed values for n < 15.
 */
double numer_ln_fact(unsigned long n) {
  if(n < 15) {
    switch(n) {
    case(0):
    case(1):
      return 0.0;
    case(2):
      return 0.6931471805599453;
    case(3):
      return 1.7917594692280550;
    case(4):
      return 3.1780538303479458;
    case(5):
      return 4.7874917427820458;
    case(6):
      return 6.579251212010101;
    case(7):
      return 8.5251613610654147;
    case(8):
      return 10.60460290274525;
    case(9):
      return 12.801827480081469;
    case(10):
      return 15.104412573075516;
    case(11):
      return 17.502307845873887;
    case(12):
      return 19.98721449566188;
    case(13):
      return 22.55216385312342;
    case(14):
      return 25.191221182738680;
    }
  }

  return (double)n * log(n) - (double)n;
}


/**
 * Returns an approximation the the factorial
 * n! = n * (n-1) * (n-2) * ... * 2 * 1
 *
 * Uses a varient of stirlings approximation for n >= 15.  
 * (see: http://mathworld.wolfram.com/StirlingsApproximation.html)
 *
 * For speed and accuracy, returns pre-computed values for n < 15.
 *
 * It is better to use the numer_ln_fact function for values of n that
 * are even somewhat large. This function will quickly exceed the
 * maximum double value.
 */ 
double numer_fact(unsigned long n) {
  if(n < 15) {
    switch(n) {
    case(0):
    case(1):
      return 1.0;
    case(2):
      return 2.0;
    case(3):
      return 6.0;
    case(4):
      return 24.0;
    case(5):
      return 120.0;
    case(6):
      return 720.0;
    case(7):
      return 5040.0;
    case(8):
      return 40320.0;
    case(9):
      return 362880.0;
    case(10):
      return 3628800.0;
    case(11):
      return 39916800.0;
    case(12):
      return 479001600.0;
    case(13):
      return 6227020800.0;
    case(14):
      return 87178291200.0;
    }
  }
  return sqrt(M_PI * (2.0*n + 1.0/3.0))* pow(n,n) * exp(-n);
}





