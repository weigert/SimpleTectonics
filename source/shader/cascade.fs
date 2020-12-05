#version 430 core
in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;
uniform sampler2D cluster;

vec3 black = vec3(0.0, 0.0, 0.0);
vec3 white = vec3(1.0, 1.0, 1.0);

uniform float D;

const float maxdiff = 0.0;

layout (std430, binding = 0) buffer segheight {
  float h[];
};


int index(vec3 a){
  a = a*255.0f;
  return int((a.x + a.y*256 + a.z*256*256));
}

vec3 cascade(vec3 p, vec3 n, float m){

  //Pile Size Difference
  vec3 diff = n - p;
  vec3 excess = (abs(diff)-m)/2.0;

//  if(diff.x < 0) excess.x = min(excess.x, n.x);
//  if(diff.y < 0) excess.y = min(excess.y, n.y);
//  if(diff.z < 0) excess.z = min(excess.z, n.z);

//  if(diff.x > 0) excess.x = min(excess.x, p.x);
//  if(diff.y > 0) excess.y = min(excess.y, p.y);
//  if(diff.z > 0) excess.z = min(excess.z, p.z);

  if(diff.x < 0) excess.x *= -1;
  if(diff.y < 0) excess.y *= -1;
  if(diff.z < 0) excess.z *= -1;

  if(diff.x == 0) excess.x = 0.0;
  if(diff.y == 0) excess.y = 0.0;
  if(diff.z == 0) excess.z = 0.0;

  return excess;
}

vec3 diffuse(){
  vec2 p = ex_Tex; //Scale
  float k = 25;

  vec3 c = textureOffset(map, p, ivec2(0,0)).xyz; //Self
  float ch = k*h[index(textureOffset(cluster, p, ivec2(0,0)).xyz)];

  vec3 a = cascade(c + ch, textureOffset(map, p, ivec2( 1, 0)).xyz + k*h[index(textureOffset(cluster, p, ivec2( 1, 0)).xyz)], maxdiff );
      a += cascade(c + ch, textureOffset(map, p, ivec2( 0, 1)).xyz + k*h[index(textureOffset(cluster, p, ivec2( 0, 1)).xyz)], maxdiff );
      a += cascade(c + ch, textureOffset(map, p, ivec2(-1, 0)).xyz + k*h[index(textureOffset(cluster, p, ivec2(-1, 0)).xyz)], maxdiff );
      a += cascade(c + ch, textureOffset(map, p, ivec2( 0,-1)).xyz + k*h[index(textureOffset(cluster, p, ivec2( 0,-1)).xyz)], maxdiff );
      //a += cascade(c + ch,textureOffset(map, p, ivec2(-1,-1)).xyz + 25*sqrt(2)*h[index(textureOffset(cluster, p, ivec2(-1,-1)).xyz)], abs(ch - h[index(textureOffset(cluster, p, ivec2(-1,-1)).xyz)]) + maxdiff );
      //a += cascade(c + ch,textureOffset(map, p, ivec2( 1,-1)).xyz + 25*sqrt(2)*h[index(textureOffset(cluster, p, ivec2( 1,-1)).xyz)], abs(ch - h[index(textureOffset(cluster, p, ivec2( 1,-1)).xyz)]) + maxdiff );
      //a += cascade(c + ch,textureOffset(map, p, ivec2(-1, 1)).xyz + 25*sqrt(2)*h[index(textureOffset(cluster, p, ivec2(-1, 1)).xyz)], abs(ch - h[index(textureOffset(cluster, p, ivec2(-1, 1)).xyz)]) + maxdiff );
      //a += cascade(c + ch,textureOffset(map, p, ivec2( 1, 1)).xyz + 25*sqrt(2)*h[index(textureOffset(cluster, p, ivec2( 1, 1)).xyz)], abs(ch - h[index(textureOffset(cluster, p, ivec2( 1, 1)).xyz)]) + maxdiff );

  return c+0.2*a;
}

void main(){
  fragColor = vec4(diffuse(), 1.0);
}
