#version 330 core

in vec2 ex_Quad;
in vec2 ex_Centroid;
flat in vec3 ex_Color;
out vec4 fragColor;

uniform float R;

float L2(vec2 a, vec2 b){
  return length(a-b);
}

void main(){
  gl_FragDepth = L2(ex_Quad, ex_Centroid);
  fragColor = vec4(ex_Color, 1.0);
}
