#version 430 core

in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;
uniform sampler2D cluster;

layout (std430, binding = 0) buffer speed {
  vec2 v[];
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

vec3 convect(float rate){

  vec2 p = ex_Tex;                                                  //Position
  int c = num(textureOffset(map, p, ivec2( 0, 0)).xyz);             //Speed
  vec2 s = v[num(textureOffset(cluster, p, ivec2( 0, 0)).xyz)];   //Height

  int cfx = num(textureOffset(map, p, ivec2( 0, 1)).xyz);
  int cbx = num(textureOffset(map, p, ivec2( 0,-1)).xyz);
  int cfy = num(textureOffset(map, p, ivec2( 1, 0)).xyz);
  int cby = num(textureOffset(map, p, ivec2(-1, 0)).xyz);

  vec2 sfx = v[num(textureOffset(cluster, p, ivec2( 0, 1)).xyz)];
  vec2 sbx = v[num(textureOffset(cluster, p, ivec2( 0,-1)).xyz)];
  vec2 sfy = v[num(textureOffset(cluster, p, ivec2( 1, 0)).xyz)];
  vec2 sby = v[num(textureOffset(cluster, p, ivec2(-1, 0)).xyz)];

  int diff = 0;                                 //Difference

  if(sfx.x > 0) diff += abs(int(sfx.x*(cfx)));  //Neighbor Contribution
  if(sbx.x < 0) diff += abs(int(sbx.x*(cbx)));
  if(sfy.y < 0) diff += abs(int(sfy.y*(cfy)));
  if(sby.y > 0) diff += abs(int(sby.y*(cby)));

  diff -= abs(int(s.x*c));                      //Self Distribution
  diff -= abs(int(s.y*c));

  return col(int(c + rate*diff));

}

void main(){
  fragColor = vec4(convect(0.015f), 1.0);
}
