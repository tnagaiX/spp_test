#version 120

uniform sampler2D QueryEmissionPhotonCountTexture;
uniform sampler2D RandomTexture;

uniform int UseEyeRays;
uniform int PathLength;
uniform vec4 BufInfo;

uniform float MaxPathLength;
uniform float Time;
uniform float Wavelength;

const bool FullSpectrum = true; // modify progressive.fs as well
const bool DOF = true;
//const bool MotionBlur = true;
const bool Glossy = true;

const float FocalLength = 210.0;
const float ApertureSize = 7.0;
const float MotionSize = 10.0;
const float Glossiness = 0.5;



// WCMYRGBXYZ using Gaussians fitting
const float White = 1.0;

const vec4 Cyan0 = vec4(0.0, 0.424537460743542, 0.0866503554583976, 0.560757618949125);
const vec4 Cyan1 = vec4(0.0, 0.246400896854156, 0.0795161416808855, 0.216116362841135);
const vec4 Cyan2 = vec4(1.0, 0.067666394964209, 0.2698588575757230, 0.890716186803857);

const vec4 Magenta0 = vec4(0.0, 0.092393363155047, -0.030670840714796, 0.425200104381996);
const vec4 Magenta1 = vec4(0.0, 0.174734179228986, 0.0690508593874629, 0.983929883263911);
const vec4 Magenta2 = vec4(1.0, 0.613995338323662, 0.0794711389383399, 1.003105061865860);

const vec4 Yellow0 = vec4(0.0, 0.369673263739623, -0.071355497310236, 0.503666150930812);
const vec4 Yellow1 = vec4(0.0, 0.558410218684172,  0.151858057162275, 0.878349029651678);
const vec4 Yellow2 = vec4(1.0, 0.587945864428471,  0.101005427723483, 0.109960421083442);

const vec4 Red0 = vec4(0.0, 0.574803873802654,  0.0349961565910619, 0.670478585641923);
const vec4 Red1 = vec4(0.0, 0.042753652345675, -0.076576978780864,  0.070884754752968);
const vec4 Red2 = vec4(1.0, 0.669048230499984,  0.0587027396330119, 0.957999219817480);

const vec4 Green0 = vec4(0.0, 0.305242141596798,  0.0337596436768638, 0.424248514020785);
const vec4 Green1 = vec4(0.0, 0.476992126451749, -0.0541085157876399, 0.815789194891182);
const vec4 Green2 = vec4(1.0, 0.365833471799225, -0.0583175076362409, 0.792406519710127);

const vec4 Blue0 = vec4(0.0, 0.144760614900738, 0.0848347582999023, 0.993361426917213);
const vec4 Blue1 = vec4(0.0, 0.600421286424602, -0.060880809655396, 0.0744873773945442);
const vec4 Blue2 = vec4(1.0, 0.231505955455338, -0.029894351908322, 0.339396172335299);



float Gaussian(const float x0, const float s, const float w, const float x)
{
  return w * exp( -(x - x0) * (x - x0) / (2.0 * s * s + 1.0e-20) );
}

float GPURnd(inout vec4 n)
{
	// from http://gpgpu.org/forums/viewtopic.php?t=2591&sid=17051481b9f78fb49fba5b98a5e0f1f3
	const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
	const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
	const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
	const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);

	vec4 beta = floor(n / q);
	vec4 p = a * (n - beta * q) - beta * r;
	beta = (sign(-p) + vec4(1.0)) * vec4(0.5) * m;

	n = (p + beta);

	return fract(dot(n / m, vec4(1.0, -1.0, 1.0, -1.0))); 
}



struct Ray
{
	vec3 org;
	vec3 dir;
};
struct Sphere
{
	vec3 c;
	float r;
	vec3 col;
	int f;
};
struct Intersection
{
	float t;
	Sphere sphere;
};



vec3 tangent(const vec3 n)
{
	vec3 t = n;
	vec3 a = abs(n);

	if ((a.x < a.y) && (a.x < a.z))
	{
		t.x = 1.0;
	}
	else if (a.y < a.z)
	{
		t.y = 1.0;
	}
	else
	{
		t.z = 1.0;
	}

	return normalize(cross(t, n));
}



vec3 along(const vec3 v, const vec3 n)
{
	vec3 t = tangent(n);
	vec3 s = cross(t, n);

	return (v.x * t + v.y * n + v.z * s);
}



vec3 glossy_normal(const vec3 n, const float g, inout vec4 rndv)
{
	float rnd1 = GPURnd(rndv);
	float rnd2 = GPURnd(rndv);

	float temp1 = 2.0 * 3.141592 * rnd1;
	float temp2 = sqrt(1.0 - pow(rnd2, 2.0 / (g + 1.0)));
	vec3 v = vec3(sin(temp1) * temp2, pow(rnd2, 1.0 / (g + 1.0)), cos(temp1) * temp2);

	return normalize(along(v, n));
}



void shpere_intersect(const Sphere sphere, const Ray ray, inout Intersection isect)
{
	vec3 rs = ray.org - sphere.c;
	float B = dot(rs, ray.dir);
	float C = dot(rs, rs) - (sphere.r * sphere.r);
	float D = B * B - C;

	if (D > 0.0)
	{
		D = sqrt(D);
		float t0 = -B - D;
		float t1 = -B + D;
		if ( (t0 > 0.0) && (t0 < isect.t) )
		{
			isect.t = t0;
			isect.sphere = sphere;
		}
		else if ( (t1 > 0.0) && (t1 < isect.t) )
		{
			isect.t = t1;
			isect.sphere = sphere;
		}
	}
}



Sphere sphere[9];
void Intersect(const Ray ray, inout Intersection isect)
{
	shpere_intersect(sphere[0], ray, isect);
	shpere_intersect(sphere[1], ray, isect);
	shpere_intersect(sphere[2], ray, isect);
	shpere_intersect(sphere[3], ray, isect);
	shpere_intersect(sphere[4], ray, isect);
	shpere_intersect(sphere[5], ray, isect);
	shpere_intersect(sphere[6], ray, isect);
	shpere_intersect(sphere[7], ray, isect);
	shpere_intersect(sphere[8], ray, isect);
}

void setScene()
{
	// scene definition
	// walls
	sphere[0].c   = vec3( 1.0e5+1.0,40.8,81.6);
	sphere[0].r   = 1.0e5;
	sphere[0].col = vec3(0.75, 0.25, 0.25);//vec3(0.75);//
	sphere[0].f = 0;

	sphere[1].c   = vec3(-1.0e5+99.0,40.8,81.6);
	sphere[1].r   = 1.0e5;
	sphere[1].col = vec3(0.25, 0.25, 0.75);//vec3(0.75);//
	sphere[1].f = 0;

	sphere[2].c   = vec3(50.0,40.8, 1.0e5);
	sphere[2].r   = 1.0e5;
	sphere[2].col = vec3(0.75);
	sphere[2].f = 0;

	sphere[3].c   = vec3(50.0,40.8,-1.0e5+170.0);
	sphere[3].r   = 1.0e5;
	sphere[3].col = vec3(0.0);
	sphere[3].f = 0;

	sphere[4].c   = vec3(50.0, 1.0e5, 81.6);
	sphere[4].r   = 1.0e5;
	sphere[4].col = vec3(0.75);
	sphere[4].f = 0;

	sphere[5].c   = vec3(50.0, -1.0e5+81.6, 81.6);
	sphere[5].r   = 1.0e5;
	sphere[5].col = vec3(0.75);
	sphere[5].f = 0;

	// spheres
	sphere[6].c   = vec3(27, 15.5, 88.0);
	sphere[6].r   = 15.5;
	sphere[6].col = vec3(0.999);
	sphere[6].f = 1;

	sphere[7].c   = vec3(73, 15.5, 88.0);
	sphere[7].r   = 15.5;
	sphere[7].col = vec3(0.999);
	sphere[7].f = 2;

	// motion blur -> no motion blur
	sphere[8].c   = vec3(50.0, 15.5, 48.0 );
	sphere[8].r   = 15.5;
	sphere[8].col = vec3(0.75);
	sphere[8].f = 0;

	vec2 PixelIndex = gl_FragCoord.xy * BufInfo.zw;

	vec4 rnd = texture2D(RandomTexture, PixelIndex);
	vec4 rndv = rnd;
}


void main()
{
    // scene definition
    setScene();

	vec2 PixelIndex = gl_FragCoord.xy * BufInfo.zw;

	vec4 rnd = texture2D(RandomTexture, PixelIndex);
	vec4 rndv = rnd;

	Ray r;
	if (!bool(UseEyeRays))
	{
        // photon trace
		// photon rays
		// point light source
		r.org = vec3(50.0, 55.0, 68.0);

		// generate random direction
		rnd.x = GPURnd(rndv);
		rnd.y = GPURnd(rndv);
		rnd.x = 2.0 * 3.141592 * rnd.x;
		rnd.y = 2.0 * acos(sqrt(1.0 - rnd.y));
		r.dir.x = sin(rnd.x) * sin(rnd.y);
		r.dir.y = cos(rnd.y);
		r.dir.z = cos(rnd.x) * sin(rnd.y);
	}
	else
	{
		// eye rays
		// pinhole camera
		const float FilmSize = 0.5135;
		const vec3 pinhole = vec3(50.0, 48.0, 295.6);
		const float AspectRatio = 1.0;

		vec3 camd = normalize(vec3(0.0, -0.042612, -1.0));

		vec3 cx = vec3(AspectRatio * FilmSize, 0.0, 0.0);
		vec3 cy = normalize(cross(cx, camd)) * FilmSize;

		vec3 d = cx * 0.5 * (PixelIndex.x * 2.0 - 1.0 + (GPURnd(rndv) - 0.5) * BufInfo.z) + 0.5 * cy * (PixelIndex.y * 2.0 - 1.0 + (GPURnd(rndv) - 0.5) * BufInfo.w) + camd;
		r.dir = normalize(d);
		r.org = pinhole + r.dir * 140.0;

		// thin lens model with a circular aperture
		if (DOF)
		{
			vec3 f = pinhole + r.dir * FocalLength;
			float radius = sqrt(GPURnd(rndv)) * ApertureSize;
			float theta = 2.0 * 3.141592 * GPURnd(rndv);
			vec3 lens = pinhole +  radius * (cx * cos(theta) + cy * sin(theta));
			r.dir = normalize(f - lens);
			r.org = lens + r.dir * 140.0;
		}
	}

	vec3 col = vec3(1.0, 1.0, 1.0);
	vec3 nrm = vec3(0.0, 0.0, 0.0);
	vec3 pos = vec3(0.0, 0.0, 0.0);
	const float eps = 0.01;

	Intersection i;
	vec3 prev_col;
	int brdf = -1;

	i.t = 1.0e+30;
	Intersect(r, i);
	if (i.t != 1.0e+30)
	{
		brdf = i.sphere.f;
		prev_col = i.sphere.col.rgb;

		pos = r.org + r.dir * i.t;
		nrm = normalize(pos - i.sphere.c);
	}

	for (int j = 1; j <= PathLength; j++)
	{
		bool total = false;

		if (brdf == 0) 
		{
			// Lambertian
			// eye ray stops at Lambertian
			if (UseEyeRays == 1)
			{
				break;
			}

			rnd.x = GPURnd(rndv);
			rnd.y = GPURnd(rndv);

			rnd.x = 2.0 * 3.141592 * rnd.x;
			rnd.y = sqrt(rnd.y);
			vec3 v = vec3(sin(rnd.x) * rnd.y, sqrt(1.0 - rnd.y * rnd.y), cos(rnd.x) * rnd.y);

			nrm = sign(-dot(nrm, r.dir)) * nrm;
			r.org = pos + eps * nrm;
			r.dir = along(v, nrm);
		}
		else if (brdf == 1)
		{
			// specular reflection
			if (Glossy)
			{
				// note that direct illumination on glossy surface is ignored in this demo
				// which should be separetely added to the final result
				nrm = glossy_normal(nrm, 1.0 / pow((1.0 - Glossiness) * 0.5, 2.71828), rndv);
			}
			nrm = sign(-dot(nrm, r.dir)) * nrm;

			r.org = pos + eps * nrm;
			r.dir = reflect(r.dir, nrm);

			if (Glossy)
			{
				if (dot(r.dir, nrm) < 0.0) r.dir = -r.dir;
			}
		}
		else if (brdf == 2)
		{
			// specular refraction
			float ln = dot(nrm, r.dir);
			float eta = 1.5220;

			float R0 = ((eta - 1.0) * (eta - 1.0)) / ((eta + 1.0) * (eta + 1.0));

			if (ln < 0.0)
			{
				// in
				float c = 1.0 + ln;
				float Re = R0 + (1.0 - R0) * c * c * c * c * c;
				float P = (Re + 1.0) / 2.0;

				if (GPURnd(rndv) < P)
				{
					prev_col *= Re / P;
					r.org = pos + eps * nrm;
					r.dir = reflect(r.dir, nrm);
				}
				else
				{
					prev_col *= (1.0 - Re) / (1.0 - P);
					r.org = pos - eps * nrm;
					r.dir = refract(r.dir, nrm, 1.0 / eta);
				}
			}
			else
			{
				// out
				float cos2t = 1.0 - eta * eta * (1.0 - ln * ln);
				if (cos2t < 0.0)
				{
					r.org = pos - eps * nrm;
					r.dir = reflect(r.dir, -nrm);
					total = true;
				}
				else
				{
					float c = 1.0 - ln;
					float Re = R0 + (1.0 - R0) * c * c * c * c * c;
					float P = (Re + 1.0) / 2.0;

					if (GPURnd(rndv) < P)
					{
						prev_col *= Re / P;
						r.org = pos - eps * nrm;
						r.dir = reflect(r.dir, -nrm);
					}
					else
					{
						prev_col *= (1.0 - Re) / (1.0 - P);
						r.org = pos + eps * nrm;
						r.dir = refract(r.dir, -nrm, eta);
					}
				}
			}
		}

		i.t = 1.0e+30;
		Intersect(r, i);
		if (i.t != 1.0e+30)
		{
			brdf = i.sphere.f;
			col *= prev_col;

			pos = r.org + r.dir * i.t;
			nrm = normalize(pos - i.sphere.c);

			prev_col = i.sphere.col.rgb;
		}
	}

	if (bool(UseEyeRays))
	{
		// eye ray tracing
        col = col * prev_col;

		if (dot(nrm, r.dir) > 0.0) nrm = -nrm;
	}
	else
	{
		// photon tracing
		col *= vec3(5000.0) * (4.0 * 3.141592) * MaxPathLength;
		nrm = r.dir;
	}

	gl_FragData[0] = vec4(pos, 1.0);
	gl_FragData[1] = vec4(col, 1.0);
	gl_FragData[2] = vec4(nrm, 1.0);
	gl_FragData[3] = rndv;
}
