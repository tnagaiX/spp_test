#version 120

uniform sampler2D QueryFluxRadiusTexture;
uniform sampler2D QueryEmissionPhotonCountTexture;

uniform float TotalPhotonNum;
uniform vec4 BufInfo;

const float gamma = 2.2;
const float exposure = 1.0;

void main()
{	
	vec2 PixelIndex = gl_FragCoord.xy * BufInfo.zw;

	vec4 QueryFluxRadius = texture2D(QueryFluxRadiusTexture, PixelIndex);

	vec3 QueryFlux = QueryFluxRadius.xyz;
	float QueryRadius = QueryFluxRadius.w;

	gl_FragColor = vec4(QueryFlux / (QueryRadius * QueryRadius * 3.141592 * TotalPhotonNum), 1.0);
	gl_FragColor = max(gl_FragColor, vec4(0.0));
	gl_FragColor = pow(vec4(1.0) - exp(-gl_FragColor * exposure), vec4(1.0 / gamma));
}
