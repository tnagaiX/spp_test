#ifndef _RANDOM_H
#define _RANDOM_H

namespace XORShift
{
	// XOR shift PRNG
	unsigned int x = 123456789;
	unsigned int y = 362436069;
	unsigned int z = 521288629;
	unsigned int w = 88675123;

	inline float frand()
	{
		const unsigned int t = x ^ (x << 11);
		x = y; y = z; z = w;
		return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8))) * (1.0f / 4294967295.0f);
	}
}

#endif