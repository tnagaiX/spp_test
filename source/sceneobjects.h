#ifndef	_SCENEOBJECTS_H_
#define	_SCENEOBJECTS_H_

#include <math.h>
#include "vec.h"
#include "ray.h"
#include "materials.h"

// シーンオブジェクトを定義する

// 球のクラス
struct Sphere {
	double rad; Vec p, c; Refl_t refl;
	Sphere(double r_, Vec p_, Vec c_, Refl_t re_) : rad(r_), p(p_), c(c_), refl(re_){}
	inline double intersect(const Ray &r) const {
		// ray-sphere intersection returns distance
		Vec op = p - r.o;
		double t, b = op.dot(r.d), det = b*b - op.dot(op) + rad*rad;
		if (det < 0) {
			return 1e20;
		}
		else {
			det = sqrt(det);
		}
		return (t = b - det) > 1e-4 ? t : ((t = b + det)>1e-4 ? t : 1e20);
	}
};


#endif