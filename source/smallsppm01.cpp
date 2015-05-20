// smallsppm01.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
#include "stdafx.h"

// expanded smallppm (code is exactly the same as smallppm.cpp but with more comments)

#include <math.h>   // smallppm, Progressive Photon Mapping by T. Hachisuka
#include <stdlib.h> // originally smallpt, a path tracer by Kevin Beason, 2008
#include <stdio.h>  // Usage: ./smallppm 100000 && xv image.ppm

#include "constant.h"
#include "vec.h"
#include "bbox.h"
#include "sceneobjects.h"
#include "random.h"
#include "scene.h"


struct HPoint {
	Vec f, pos, nrm, flux;
	double r2;
	unsigned int n; // n = N / ALPHA in the paper
	int pix;
};

struct List { HPoint *id; List *next; };
List* ListAdd(HPoint *i, List* h){
	List* p = new List;
	p->id = i;
	p->next = h;
	return p;
}

unsigned int num_hash, pixel_index, num_photon;
double hash_s; List **hash_grid; List *hitpoints = NULL; AABB hpbbox;

// spatial hash function
inline unsigned int hash(const int ix, const int iy, const int iz) {
	return (unsigned int)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % num_hash;
}

void build_hash_grid(const int w, const int h) {
	// find the bounding box of all the measurement points
	hpbbox.reset();
	List *lst = hitpoints;
	while (lst != NULL) {
		HPoint *hp = lst->id;
		lst = lst->next;
		hpbbox.fit(hp->pos);
	}

	// heuristic for initial radius
	Vec ssize = hpbbox.max - hpbbox.min;
	double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;

	// determine hash table size
	// we now find the bounding box of all the measurement points inflated by the initial radius
	hpbbox.reset();
	lst = hitpoints;
	int vphoton = 0;
	while (lst != NULL) {
		HPoint *hp = lst->id;
		lst = lst->next;
		hp->r2 = irad *irad;
		hp->n = 0;
		hp->flux = Vec();
		vphoton++;
		hpbbox.fit(hp->pos - irad);
		hpbbox.fit(hp->pos + irad);
	}

	// make each grid cell two times larger than the initial radius
	hash_s = 1.0 / (irad*2.0);
	num_hash = vphoton;

	// build the hash table
	hash_grid = new List*[num_hash];
	for (unsigned int i = 0; i<num_hash; i++) hash_grid[i] = NULL;
	lst = hitpoints;
	while (lst != NULL)
	{
		HPoint *hp = lst->id;
		lst = lst->next;
		Vec BMin = ((hp->pos - irad) - hpbbox.min) * hash_s;
		Vec BMax = ((hp->pos + irad) - hpbbox.min) * hash_s;
		for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
		{
			for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
			{
				for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
				{
					int hv = hash(ix, iy, iz);
					hash_grid[hv] = ListAdd(hp, hash_grid[hv]);
				}
			}
		}
	}
}

// tone mapping and gamma correction
int toInt(double x){
	return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5);
}

// find the closet interection
inline bool intersect(const Ray &r, double &t, int &id){
	int n = sizeof(sph) / sizeof(Sphere);
	double d, inf = 1e20; t = inf;
	for (int i = 0; i<n; i++){
		d = sph[i].intersect(r);
		if (d<t){
			t = d;
			id = i;
		}
	}
	return t<inf;
}

// generate a photon ray from the point light source with QMC
// フォトンの生成
void generate_photon(Ray* ray_ptr, Vec* f, int i) {
	// フォトン１個あたりのフラックス、変更が必要？
	// 4πで割っていて、角度について乱数で振っているので
	// 立体角測度に従って、モンテカルロ積分をしているのに相当する
	*f = Vec(2500, 2500, 2500)*(PI*4.0); 

	// 確率的に、φとθを計算
	// これ、完全拡散面のインポータンスサンプリングのような形になっているが、これでよいのか？
	double p = 2.*PI*hal(0, i);
	double t = 2.*acos(sqrt(1. - hal(1, i)));

	double st = sin(t);

	// 方向のばらつきは、φとθのみとなっているので
	// このままでも使えるのでは？
	// oをそのまま飛ばすと、cornel_boxの外に行ってしまうので
	// 画角を制限するほうがよいか？
	//ray_ptr->d = Vec(cos(p)*st, cos(t), sin(p)*st);
	ray_ptr->o = Vec(50, 60, 85);

	double phi = 2.0 * PI * hal(0, i);
	double z = hal(1, i) * 2.0 - 1.0;

	double sz = sqrt(1.0 - z*z);

	double aaa = hal(1, i);

	ray_ptr->d = Vec(sz * cos(phi), sz * sin(phi), z);
}

void trace(const Ray &r, int dpt, bool m, const Vec &fl, const Vec &adj, int i)
{
	double t;
	int id;

	dpt++;
	if (!intersect(r, t, id) || (dpt >= 20))return;

	int d3 = dpt * 3;
	const Sphere &obj = sph[id];
	Vec x = r.o + r.d*t, n = (x - obj.p).norm(), f = obj.c;
	Vec nl = n.dot(r.d)<0 ? n : n*-1;
	double p = f.x>f.y&&f.x>f.z ? f.x : f.y>f.z ? f.y : f.z;

	if (obj.refl == DIFF) {
		// Lambertian

		// use QMC to sample the next direction
		double r1 = 2.*PI*hal(d3 - 1, i), r2 = hal(d3 + 0, i);
		double r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x)>.1 ? Vec(0, 1) : Vec(1)).cross(w)).norm();
		Vec v = w.cross(u), d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();

		if (m) {
			// eye ray
			// store the measurment point
			HPoint* hp = new HPoint;
			hp->f = f.mul(adj);
			hp->pos = x;
			hp->nrm = n;
			hp->pix = pixel_index;
			hitpoints = ListAdd(hp, hitpoints);
		}
		else
		{
			// photon ray
			// find neighboring measurement points and accumulate flux via progressive density estimation
			Vec hh = (x - hpbbox.min) * hash_s;
			int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
			// strictly speaking, we should use #pragma omp critical here.
			// it usually works without an artifact due to the fact that photons are 
			// rarely accumulated to the same measurement points at the same time (especially with QMC).
			// it is also significantly faster.
			{
				List* hp = hash_grid[hash(ix, iy, iz)];
				while (hp != NULL) {
					HPoint *hitpoint = hp->id;
					hp = hp->next;
					Vec v = hitpoint->pos - x;
					// check normals to be closer than 90 degree (avoids some edge brightning)
					if ((hitpoint->nrm.dot(n) > 1e-3) && (v.dot(v) <= hitpoint->r2)) {
						// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
						double g = (hitpoint->n*ALPHA + ALPHA) / (hitpoint->n*ALPHA + 1.0);
						hitpoint->r2 = hitpoint->r2*g;
						hitpoint->n++;
						hitpoint->flux = (hitpoint->flux + hitpoint->f.mul(fl)*(1. / PI))*g;
					}
				}
			}
			if (hal(d3 + 1, i)<p) trace(Ray(x, d), dpt, m, f.mul(fl)*(1. / p), adj, i);
		}

	}
	else if (obj.refl == SPEC) {
		// mirror
		trace(Ray(x, r.d - n*2.0*n.dot(r.d)), dpt, m, f.mul(fl), f.mul(adj), i);

	}
	else {
		// glass
		Ray lr(x, r.d - n*2.0*n.dot(r.d));
		bool into = (n.dot(nl)>0.0);
		double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

		// total internal reflection
		if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0) return trace(lr, dpt, m, fl, adj, i);

		Vec td = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
		double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : td.dot(n));
		double Re = R0 + (1 - R0)*c*c*c*c*c, P = Re; Ray rr(x, td); Vec fa = f.mul(adj);
		if (m) {
			// eye ray (trace both rays)
			trace(lr, dpt, m, fl, fa*Re, i);
			trace(rr, dpt, m, fl, fa*(1.0 - Re), i);
		}
		else {
			// photon ray (pick one via Russian roulette)
			(hal(d3 - 1, i)<P) ? trace(lr, dpt, m, fl, fa, i) : trace(rr, dpt, m, fl, fa, i);
		}
	}
}

// eyeレイトレーシング（名称はいずれ見直すかも）
void eye_ray_tracing(int width, int height){
	// trace eye rays and store measurement points
	//Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
	Ray cam(Vec(50, 48, 220.0), Vec(0, -0.042612, -1).norm());

	double focal_length = 140.0;

	// イメージセンサの縦幅（物理的大きさ）
	double imagesensor_size = 30.0;
	double sensor_to_lens_dist = 30.0;

	double imagesensor_width = imagesensor_size * (double)width / (double)height;
	double imagesensor_height = imagesensor_size;

	// イメージセンサの上方向
	Vec imagesensor_up = Vec(0.0, 1.0, 0.0);

	// イメージセンサのuvを求める
	Vec imagesensor_u, imagesensor_v;
	imagesensor_u = cam.d.cross(imagesensor_up).norm() * imagesensor_width;
	imagesensor_v = imagesensor_u.cross(cam.d) * imagesensor_height;

	// オブジェクトプレーンの計算
	Vec objplane_center = cam.o + cam.d * (focal_length + sensor_to_lens_dist);
	Vec objplane_u = imagesensor_u;
	Vec objplane_v = imagesensor_v;

	double sensor_size = 0.5135;

	// コードの書き方として、いったんcamの位置にイメージセンサが置かれている、として書く
	// これがeduptのimagesensor_sensor相当？
	// これに合わせて、scene.hも書き直し

	Vec cx = Vec(width * sensor_size / height);

	// 外積の演算は書き換えたい
	Vec cy = (cx.cross(cam.d)).norm()*sensor_size;


	// eyeパスのトレーシング
	for (int y = 0; y < height; y++){
		fprintf(stderr, "\rEyePass %5.2f%%", 100.0*y / (height - 1));
		for (int x = 0; x < width; x++) {
			// イメージセンサ上の位置を計算
			// [-0.5, 0.5] (0,0)が中心
			double u_on_imagesensor = (double)x / (double)width - 0.5;
			double v_on_imagesensor = (double)y / (double)height - 0.5;

			Vec pos_on_imagesensor = cam.o
				+ imagesensor_u * u_on_imagesensor
				+ imagesensor_v * v_on_imagesensor;

			// オブジェクトプレーン上の

			//pixel_index = x + y * width;
			//Vec d = cx * ((x + 0.5) / width - 0.5) + cy * (-(y + 0.5) / height + 0.5) + cam.d;



			trace(Ray(cam.o + d * 140, d.norm()), 0, true, Vec(), Vec(1, 1, 1), 0);
		}
	}
	fprintf(stderr, "\n");

	// build the hash table over the measurement points
	build_hash_grid(width, height);
}



int main(int argc, char *argv[]) {
	// samps * 1000 photon paths will be traced
	int width = 640;
	int height = 480;
	int samples = (argc == 2) ? MAX(atoi(argv[1]) / 1000, 1) : 1000;

	// eyeパスのトレーシング
	eye_ray_tracing(width, height);
	
	// trace photon rays with multi-threading
	num_photon = samples;

	// バッファの宣言
	Vec *c = new Vec[width * height];
	Vec vw = Vec(1, 1, 1);

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i<num_photon; i++) {
		double p = 100.*(i + 1) / num_photon;
		fprintf(stderr, "\rPhotonPass %5.2f%%", p);
		int m = 1000 * i;
		Ray ray;
		Vec flux;

		// To Do縮減の繰り返し回数を指定できるようにする
		for (int j = 0; j<1000; j++){
			generate_photon(&ray, &flux, m + j);
			trace(ray, 0, 0>1, flux, vw, m + j);
		}
	}

	// density estimation
	List* lst = hitpoints;
	while (lst != NULL) {
		HPoint* hp = lst->id;
		lst = lst->next;
		int i = hp->pix;
		c[i] = c[i] + hp->flux*(1.0 / (PI*hp->r2*num_photon*1000.0));
	}

	// save the image after tone mapping and gamma correction
	FILE* f = fopen("image.ppm", "w"); fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
	for (int i = 0; i< width * height; i++) {
		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	}
}