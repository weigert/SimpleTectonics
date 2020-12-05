#version 430 core
in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;
uniform sampler2D cluster;

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

vec3 convect(){

  //Position, Height, Speed
  vec2 p = ex_Tex;

  int c = num(textureOffset(map, p, ivec2( 0, 0)).xyz);             //Self
  vec2 s = v[index(textureOffset(cluster, p, ivec2( 0, 0)).xyz)];

  int cfx = num(textureOffset(map, p, ivec2( 0, 1)).xyz); //Self
  int cbx = num(textureOffset(map, p, ivec2( 0,-1)).xyz); //Self
  int cfy = num(textureOffset(map, p, ivec2( 1, 0)).xyz); //Self
  int cby = num(textureOffset(map, p, ivec2(-1, 0)).xyz); //Self

  vec2 sfx = (v[index(textureOffset(cluster, p, ivec2( 0, 1)).xyz)]); //Self
  vec2 sbx = (v[index(textureOffset(cluster, p, ivec2( 0,-1)).xyz)]); //Self
  vec2 sfy = (v[index(textureOffset(cluster, p, ivec2( 1, 0)).xyz)]); //Self
  vec2 sby = (v[index(textureOffset(cluster, p, ivec2(-1, 0)).xyz)]); //Self

  //Difference
  int diff = 0;

  //Flux that the cell gives away
  //diff -= int(c*length(s));

  //Flux that the neighbors give!
/*  diff -= int(0.5*sfx.x*cfx);
  diff -= int(0.5*sbx.x*cbx);
  diff -= int(0.5*sfy.y*cfy);
  diff -= int(0.5*sby.y*cby);
  */
  if(sfy.x < 0) diff += abs(int(sfx.x*cfx));
  if(sby.x > 0) diff += abs(int(sbx.x*cbx));
  if(sfy.y < 0) diff += abs(int(sfy.y*cfy));
  if(sby.y > 0) diff += abs(int(sby.y*cby));

  diff -= abs(int(s.x*c));
  diff -= abs(int(s.y*c));

  return col(int(c + 0.08*diff));

/*
  //Position, Height, Speed
  vec2 p = ex_Tex;

  int c = num(textureOffset(map, p, ivec2( 0, 0)).xyz);             //Self
  vec2 s = v[index(textureOffset(cluster, p, ivec2( 0, 0)).xyz)];

  int cfx = num(textureOffset(map, p, ivec2( 0, 1)).xyz); //Self
  int cbx = num(textureOffset(map, p, ivec2( 0,-1)).xyz); //Self
  int cfy = num(textureOffset(map, p, ivec2( 1, 0)).xyz); //Self
  int cby = num(textureOffset(map, p, ivec2(-1, 0)).xyz); //Self

  vec2 sfx = v[index(textureOffset(cluster, p, ivec2( 0, 1)).xyz)]; //Self
  vec2 sbx = v[index(textureOffset(cluster, p, ivec2( 0,-1)).xyz)]; //Self
  vec2 sfy = v[index(textureOffset(cluster, p, ivec2( 1, 0)).xyz)]; //Self
  vec2 sby = v[index(textureOffset(cluster, p, ivec2(-1, 0)).xyz)]; //Self

  //Difference
  int diff = 0;

  //Flux that the cell gives away
  diff -= int(c*length(s));

  //Flux that the neighbors give!
  if(sfx.x < 0) diff += abs(int(sfx.x*cfx));
  if(sbx.x > 0) diff += abs(int(sbx.x*cbx));
  if(sfy.y < 0) diff += abs(int(sfy.y*cfy));
  if(sby.y > 0) diff += abs(int(sby.y*cby));

  return col(int(c + 0.001*diff));
  */
}

void main(){
  fragColor = vec4(convect(), 1.0);
}
