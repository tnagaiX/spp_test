#version 120

uniform sampler2D PhotonIndexTexture;

uniform float PhotonBufferSize;
uniform vec4 BufInfo;

varying vec4 p;

void main()
{
	// read hashed photon index
	vec2 TexCoord = (gl_Vertex.xy + vec2(1.0)) * 0.5;
	vec4 PhotonIndex = texture2D(PhotonIndexTexture, TexCoord);
	vec2 PhotonListIndex = PhotonIndex.zw; 

	gl_Position = vec4(PhotonListIndex * BufInfo.zw * 2.0 - vec2(1.0), 0.5, 1.0);

	p = PhotonIndex;
}
