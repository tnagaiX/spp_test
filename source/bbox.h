#ifndef	_BBOX_H_
#define	_BBOX_H_

#include "vec.h"

#define MAX(x, y) ((x > y) ? x : y)

// Axis Aligned Bounding Box 
struct AABB {
	Vec min, max; // axis aligned bounding box
	inline void fit(const Vec &p)
	{
		if (p.x < min.x)min.x = p.x; // min
		if (p.y < min.y)min.y = p.y; // min
		if (p.z < min.z)min.z = p.z; // min
		max.x = MAX(p.x, max.x);
		max.y = MAX(p.y, max.y);
		max.z = MAX(p.z, max.z);
	}
	inline void reset() {
		min = Vec(1e20, 1e20, 1e20);
		max = Vec(-1e20, -1e20, -1e20);
	}
}; 


#endif