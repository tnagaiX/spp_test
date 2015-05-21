#ifndef	_SCENE_H_
#define	_SCENE_H_

#include "vec.h"
#include "materials.h"
#include "sceneobjects.h"

Sphere sph[] = { // Scene: radius, position, color, material
	Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25, .25), DIFF),//Left
	Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(.25, .25, .75), DIFF),//Right
	Sphere(1e5, Vec(50, 40.8, 1e5), Vec(.75, .75, .75), DIFF),//Back
	//Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), DIFF),//Front
	Sphere(1e5, Vec(50, 40.8, -1e5 + 250), Vec(), DIFF),//Front
	Sphere(1e5, Vec(50, 1e5, 81.6), Vec(.75, .75, .75), DIFF),//Bottomm
	Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(.75, .75, .75), DIFF),//Top
	Sphere(16.5, Vec(27, 16.5, 47), Vec(1, 1, 1)*.999, SPEC),//Mirror
	Sphere(16.5, Vec(73, 16.5, 88), Vec(1, 1, 1)*.999, REFR),//Glass
	Sphere(8.5, Vec(50, 8.5, 60), Vec(1, 1, 1)*.999, DIFF) //Middle
};



#endif