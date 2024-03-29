#version 430 core

in vec2 ex_Quad;
in vec2 ex_Centroid;
flat in int ex_InstanceID;
flat in vec3 ex_Color;
out vec4 fragColor;

uniform float R;
uniform bool depthmap;

layout (std430, binding = 0) buffer area {
  int a[];
};

float L2(vec2 a, vec2 b){
  return length(a-b);
}

void main(){
  gl_FragDepth = L2(ex_Quad, ex_Centroid);
  if(gl_FragDepth > R) discard;

  if(depthmap) fragColor = vec4(vec3(gl_FragDepth/R), 1.0);
  else  fragColor = vec4(ex_Color, 1.0);

  atomicAdd(a[ex_InstanceID], 1); //Accumulate Area
}
