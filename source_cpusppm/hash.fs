#version 120

uniform sampler2D PhotonPositionTexture;

uniform float HashNum;
uniform float HashScale;
uniform vec3 BBoxMin;

uniform vec4 BufInfo;	



vec2 convert1Dto2D(const float t)
{
	vec2 tmp;
	tmp.x = mod(t, BufInfo.x) + 0.5;
	tmp.y = floor(t * BufInfo.z) + 0.5;

	return tmp;
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


void main()
{
	vec2 PixelIndex = gl_FragCoord.xy * BufInfo.zw;

	vec4 PhotonIndex = vec4(PixelIndex.xy, -1.0, 0.0);
	vec3 PhotonPosition = texture2D(PhotonPositionTexture, PhotonIndex.xy).xyz;

	vec3 HashIndex = floor((PhotonPosition - BBoxMin) * HashScale);
	PhotonIndex.zw = convert1Dto2D(hash(HashIndex));

	gl_FragColor = PhotonIndex;
}


