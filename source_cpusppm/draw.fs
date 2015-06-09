#version 120

uniform sampler2D input_tex;
uniform vec4 BufInfo;

void main()
{
	gl_FragColor = texture2D(input_tex, gl_FragCoord.st * BufInfo.zw);
}
