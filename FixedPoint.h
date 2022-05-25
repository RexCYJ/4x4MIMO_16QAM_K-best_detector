// define a fixed-point type structure

#include <stdint.h>

#define FRAC 13
#define INTE 2
#define WLEN 15

typedef struct FixedPoint {
	int16_t val;
	FXP operator+(FXP &);
	FXP operator-(FXP &);
	FXP operator*(FXP &);
} FXP;

FXP FXP::operator+(FXP &other) {
	int32_t ans;
	ans = (int32_t)val + (int32_t)other.val;
	if ((ans >> WLEN) ^ (ans >> (WLEN + 1)))
		if (ans >> WLEN)	ans = (1 << WLEN) - 1;
		else	ans = (1 << (WLEN + 1)) - 1;
}

