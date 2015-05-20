#ifndef	_RAY_H_
#define	_RAY_H_

#include "vec.h"

struct Ray {
	Vec o;
	Vec d;

	Ray(){};
	Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

#endif