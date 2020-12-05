#version 430 core
in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;

vec3 black = vec3(0.0, 0.0, 0.0);
vec3 white = vec3(1.0, 1.0, 1.0);

const float maxdiff = 0.0;

layout (std430, binding = 0) buffer speed {
  vec2 v[];
};

int index(vec3 a){
  a = a*255.0f;
  return int((a.x + a.y*256 + a.z*256*256));
}

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

int cascade(int p, int n){

  int diff = n - p;
  int excess = int(abs(diff)/2.0f)-250000;

  if(diff == 0) return 0;
  if(diff > 0) return min(excess,n);
  if(diff < 0) return -min(excess,p);

}

vec3 diffuse(){

  vec2 p = ex_Tex; //Scale

  //Get the value at Position
  int c = num(textureOffset(map, p, ivec2(0,0)).xyz);

  int a  = cascade(c, num(textureOffset(map, p, ivec2( 1, 0)).xyz));
      a += cascade(c, num(textureOffset(map, p, ivec2(-1, 0)).xyz));
      a += cascade(c, num(textureOffset(map, p, ivec2( 0, 1)).xyz));
      a += cascade(c, num(textureOffset(map, p, ivec2( 0,-1)).xyz));
      a += cascade(c, num(textureOffset(map, p, ivec2(-1,-1)).xyz));
      a += cascade(c, num(textureOffset(map, p, ivec2( 1,-1)).xyz));
      a += cascade(c, num(textureOffset(map, p, ivec2(-1, 1)).xyz));
      a += cascade(c, num(textureOffset(map, p, ivec2( 1, 1)).xyz));

  //Return the Colorized Version
  return col(int(c+0.05*a));
}

void main(){
  fragColor = vec4(diffuse(), 1.0);
}
