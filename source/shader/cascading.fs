#version 430 core

in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;
uniform sampler2D cluster;

layout (std430, binding = 0) buffer height {
  float h[];
};

int num(vec3 a){
  a = a*255.0f;
  return int((a.x + a.y*256 + a.z*256*256));
}

vec3 col(int i){
  float r = ((i >>  0) & 0xff)/255.0f;
  float g = ((i >>  8) & 0xff)/255.0f;
  float b = ((i >> 16) & 0xff)/255.0f;
  return vec3(r,g,b);
}

int pile(int p, int n){

  const int maxdiff = 25000;

  int diff = n - p;
  int excess = int(abs(diff)/2.0f)-maxdiff;

  if(abs(diff) < maxdiff) return diff/2;
  if(diff > 0) return min(excess,n);
  if(diff < 0) return -min(excess,p);

}

vec3 cascade(float rate){

  int c = num(textureOffset(map, ex_Tex, ivec2(0,0)).xyz);

  int a  = pile(c, num(textureOffset(map, ex_Tex, ivec2( 1, 0)).xyz));
      a += pile(c, num(textureOffset(map, ex_Tex, ivec2(-1, 0)).xyz));
      a += pile(c, num(textureOffset(map, ex_Tex, ivec2( 0, 1)).xyz));
      a += pile(c, num(textureOffset(map, ex_Tex, ivec2( 0,-1)).xyz));

  return col(int(c+rate*a));

}

void main(){
  fragColor = vec4(cascade(0.25f), 1.0);
}
