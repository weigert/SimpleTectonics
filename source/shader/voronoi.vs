#version 330 core
layout (location = 0) in vec2 in_Quad;
layout (location = 2) in vec2 in_Centroid;

out vec2 ex_Quad;
flat out vec3 ex_Color;
flat out int ex_InstanceID;

out vec2 ex_Centroid;

uniform float R;

vec3 color(int i){
  float r = ((i >>  0) & 0xff);
  float g = ((i >>  8) & 0xff);
  float b = ((i >> 16) & 0xff);
  return vec3(r,g,b)/255.0f;
}

void main(){
  ex_Centroid = in_Centroid/128.0f-1.0f;
  ex_Color = color(gl_InstanceID);
  ex_InstanceID = gl_InstanceID;
  ex_Quad =  R*in_Quad+ex_Centroid;
  gl_Position = vec4(ex_Quad, 0.0, 1.0);
}
