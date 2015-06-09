#version 120

uniform sampler2D QueryPositionTexture;
uniform sampler2D QueryEmissionPhotonCountTexture;
uniform sampler2D QueryFluxRadiusTexture;
uniform sampler2D QueryReflectanceTexture;
uniform sampler2D QueryDirectionTexture;

uniform sampler2D HashedPhotonTexture;
uniform sampler2D PhotonFluxTexture;
uniform sampler2D PhotonPositionTexture;
uniform sampler2D PhotonDirectionTexture;
uniform sampler2D PhotonCorrectionTexture;

uniform vec4 BufInfo;

uniform float HashNum;

uniform float HashScale;
uniform vec3 BBoxMin;

vec3 Flux = vec3(0.0);
float PhotonCount = 0.0;
uniform float Wavelength;
const bool FullSpectrum = false; // modify raytrace.fs as well



// CIE response curves using Gaussians fitting
const vec4 CIEX0 = vec4(0.26714125, 0.173056848724526, -0.0517890668554628, 0.369341509681465);
const vec4 CIEX1 = vec4(0.0, 0.510852785928701, 0.636521548441552, -0.324530476950362);
const vec4 CIEX2 = vec4(1.0622, 0.547302197035226, 0.0899535691555178, 1.10399973088081);

const vec4 CIEY0 = vec4(0.2671425, 0.86798560108836, 0.150307921271593, -0.354744089805774);
const vec4 CIEY1 = vec4(0, 0.10539332389757, 0.168752691961971, -0.289650515359526);
const vec4 CIEY2 = vec4(1.0002, 0.445956775505726, 0.0920541376951253, 0.814888040084084);

const vec4 CIEZ0 = vec4(0.26714375, 0.174251742295476, -0.0569218355789753, 1.72408897831517);
const vec4 CIEZ1 = vec4(0.0, 0.0542544622978704, 0.0457454482464726, -0.442679263574661);
const vec4 CIEZ2 = vec4(1.7826, 0.711309229610584, 0.285040831286585, -0.407629686738774);



float Gaussian(const float x0, const float s, const float w, const float x)
{
  return w * exp( -(x - x0) * (x - x0) / (2.0 * s * s + 1.0e-20) );
}



float GaussianMixture(const float lambda, const vec4 Data0, const vec4 Data1, const vec4 Data2)
{
	float t = (lambda - 0.380) / (0.780 - 0.380);
	float g0 = Gaussian(Data0.y, Data0.z, Data0.w, t);
	float g1 = Gaussian(Data1.y, Data1.z, Data1.w, t);
	float g2 = Gaussian(Data2.y, Data2.z, Data2.w, t);

	return min(max(g0 + g1 + g2 + Data0.x, Data1.x), Data2.x);
}




vec3 Spectrum2RGB(const float lambda)
{
	float x = GaussianMixture(lambda, CIEX0, CIEX1, CIEX2);
	float y = GaussianMixture(lambda, CIEY0, CIEY1, CIEY2);
	float z = GaussianMixture(lambda, CIEZ0, CIEZ1, CIEZ2);

  // E to D65
  // 0.26713798 is for mapping spectrum 1.0 into rgb (1.0, 1.0, 1.0)
  x = x * 0.9504700 / 0.26713798;
  y = y * 1.0000000 / 0.26713798;
  z = z * 1.0888300 / 0.26713798;

  // sRGB (D65)
	vec3 rgb;
  rgb.r = (x * ( 3.2404542) + y * (-1.5371385) + z * (-0.4985314));
  rgb.g = (x * (-0.9692660) + y * ( 1.8760108) + z * ( 0.0415560));
  rgb.b = (x * ( 0.0556434) + y * (-0.2040259) + z * ( 1.0572252));

	return rgb;
}



float hash(const vec3 idx)
{
	// use the same procedure as GPURnd
	vec4 n = vec4(idx, idx.x + idx.y - idx.z) * 4194304.0 / HashScale;

	const vec4 q = vec4(   1225.0,    1585.0,    2457.0,    2098.0);
	const vec4 r = vec4(   1112.0,     367.0,      92.0,     265.0);
	const vec4 a = vec4(   3423.0,    2646.0,    1707.0,    1999.0);
	const vec4 m = vec4(4194287.0, 4194277.0, 4194191.0, 4194167.0);

	vec4 beta = floor(n / q);
	vec4 p = a * (n - beta * q) - beta * r;
	beta = (sign(-p) + vec4(1.0)) * vec4(0.5) * m;
	n = (p + beta);

	return floor( fract(dot(n / m, vec4(1.0, -1.0, 1.0, -1.0))) * HashNum );
}



vec2 convert1Dto2D(const float t)
{
	vec2 tmp;
	tmp.x = mod(t, BufInfo.x) + 0.5;
	tmp.y = floor(t * BufInfo.z) + 0.5;

	return tmp;
}



void AccumulatePhotons(const vec3 QueryPosition, const vec3 QueryDirection, const float QueryRadius, const vec3 HashIndex)
{
	// get the first photon
	vec2 HashedPhotonIndex = convert1Dto2D(hash(HashIndex)) * BufInfo.zw;
	vec2 PhotonIndex = texture2D(HashedPhotonTexture, HashedPhotonIndex).xy;

	// accumulate photon
	vec3 PhotonPosition = texture2D(PhotonPositionTexture, PhotonIndex).xyz;
	vec3 PhotonDirection = texture2D(PhotonDirectionTexture, PhotonIndex).xyz;

	// make sure that the photon is actually in the given grid cell
	vec3 RangeMin = HashIndex / HashScale + BBoxMin;
	vec3 RangeMax = (HashIndex + vec3(1.0)) / HashScale + BBoxMin;
	if ((RangeMin.x < PhotonPosition.x) && (PhotonPosition.x < RangeMax.x) &&
	(RangeMin.y < PhotonPosition.y) && (PhotonPosition.y < RangeMax.y) &&
	(RangeMin.z < PhotonPosition.z) && (PhotonPosition.z < RangeMax.z))
	{
		float d = distance(PhotonPosition, QueryPosition);
		if ((d < QueryRadius) && (-dot(QueryDirection, PhotonDirection) > 0.001))
		{
			vec3 PhotonFlux = texture2D(PhotonFluxTexture, PhotonIndex).rgb;
			float PhotonCorrection = texture2D(PhotonCorrectionTexture, HashedPhotonIndex).x;

			Flux += PhotonFlux * PhotonCorrection;
			PhotonCount += PhotonCorrection;
		}
	}
}



void main()
{	
	vec2 PixelIndex = gl_FragCoord.xy * BufInfo.zw;

	vec3 QueryPosition = texture2D(QueryPositionTexture, PixelIndex).xyz;
	vec3 QueryDirection = texture2D(QueryDirectionTexture, PixelIndex).xyz;

	vec4 QueryFluxRadius = texture2D(QueryFluxRadiusTexture, PixelIndex);
	vec4 QueryEmissionPhotonCount = texture2D(QueryEmissionPhotonCountTexture, PixelIndex);

	vec3 QueryFlux = QueryFluxRadius.xyz;
	vec3 QueryReflectance = texture2D(QueryReflectanceTexture, PixelIndex).xyz;

	float QueryPhotonCount = QueryEmissionPhotonCount.w;
	float QueryRadius = QueryFluxRadius.w;

	vec3 RangeMin = abs(QueryPosition - vec3(QueryRadius) - BBoxMin) * HashScale;
	vec3 RangeMax = abs(QueryPosition + vec3(QueryRadius) - BBoxMin) * HashScale;

	for (int iz = int(RangeMin.z); iz <= int(RangeMax.z); iz ++)
	{
		for (int iy = int(RangeMin.y); iy <= int(RangeMax.y); iy++)
		{
			for (int ix = int(RangeMin.x); ix <= int(RangeMax.x); ix++)
			{
				AccumulatePhotons(QueryPosition, QueryDirection, QueryRadius, vec3(ix, iy, iz));
			}
		}
	}

	// BRDF (Lambertian)
	Flux = Flux * (QueryReflectance / 3.141592);
	if (FullSpectrum) Flux = Spectrum2RGB(Wavelength) * Flux.r;

	// progressive refinement
	const float alpha = 0.8;
	float g = min((QueryPhotonCount + alpha * PhotonCount) / (QueryPhotonCount + PhotonCount), 1.0);
	QueryRadius = QueryRadius * sqrt(g);
	QueryPhotonCount = QueryPhotonCount + PhotonCount * alpha;
	QueryFlux = (QueryFlux + Flux) * g;

	// output
	gl_FragData[0] = vec4(QueryFlux, QueryRadius);
	gl_FragData[1] = vec4(QueryEmissionPhotonCount.rgb, QueryPhotonCount);
}


