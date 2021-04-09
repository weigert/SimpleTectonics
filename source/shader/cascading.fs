#version 430 core

in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;
uniform sampler2D cluster;

uniform bool init;

layout (std430, binding = 0) buffer height {
  float h[];
};

int num(vec3 a){
  a = a*255.0f;
  return int((int(a.x) << 0) + (int(a.y) << 8) + (int(a.z) << 16));
}

vec3 col(int i){
  float r = ((i >>  0) & 0xff);
  float g = ((i >>  8) & 0xff);
  float b = ((i >> 16) & 0xff);
  return vec3(r,g,b)/255.0f;
}

int pile(int p, int n){

  const int maxdiff = 10000;

  int diff = n - p;
  if(abs(diff) < maxdiff) return diff/3;

  int excess = int(abs(diff)/2.0f)-maxdiff/2;

  if(diff > 0) return min(excess,n);
  if(diff < 0) return -min(excess,p);

}

vec3 cascade(float rate){

vec3 color = textureOffset(cluster, ex_Tex, ivec2(0,0)).xyz;
int c = int(255.0f*255.0f*255.0f*h[num(color)]);

if(color == vec3(1)) return col(c);

vec3 color0 = textureOffset(cluster, ex_Tex, ivec2( 1, 0)).xyz;
vec3 color1 = textureOffset(cluster, ex_Tex, ivec2(-1, 0)).xyz;
vec3 color2 = textureOffset(cluster, ex_Tex, ivec2( 0, 1)).xyz;
vec3 color3 = textureOffset(cluster, ex_Tex, ivec2( 0,-1)).xyz;

int c0 = (color0 != vec3(1))?int(255.0f*255.0f*255.0f*h[num(color0)]):0;
int c1 = (color1 != vec3(1))?int(255.0f*255.0f*255.0f*h[num(color1)]):0;
int c2 = (color2 != vec3(1))?int(255.0f*255.0f*255.0f*h[num(color2)]):0;
int c3 = (color3 != vec3(1))?int(255.0f*255.0f*255.0f*h[num(color3)]):0;

if(c < 0) c = 0;
if(c0 < 0) c0 = 0;
if(c1 < 0) c1 = 0;
if(c2 < 0) c2 = 0;
if(c3 < 0) c3 = 0;

int ah = pile(c, c0) + pile(c, c1) + pile(c, c2) + pile(c, c3);

return col(int(c+rate*ah));

}

vec3 heightcol(ivec2 offset){

  vec3 color = textureOffset(cluster, ex_Tex, ivec2(0,0)).xyz;
  if(color == vec3(1)) return vec3(0);
  int c = int(255.0f*255.0f*255.0f*h[num(color)]);
  if(c < 0.0) return col(int(0));
  return col(int(c));

}

void main(){
  //fragColor = vec4(heightcol(ivec2(0,0)), 1.0);
  if(init) fragColor = vec4(heightcol(ivec2(0,0)), 1.0);
  else fragColor = vec4(cascade(0.5f), 1.0);
}
