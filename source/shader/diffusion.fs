#version 330 core
in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;

vec3 black = vec3(0.0, 0.0, 0.0);
vec3 white = vec3(1.0, 1.0, 1.0);

float D = 0.2;

vec3 diffuse(){
  vec2 p = ex_Tex; //Scale
  vec3 c = textureOffset(map, p, ivec2(0,0)).xyz; //Self

  vec3 fx = textureOffset(map, p, ivec2( 1, 0)).xyz-c;
      fx += textureOffset(map, p, ivec2(-1, 0)).xyz-c;

  vec3 fy = textureOffset(map, p, ivec2( 0, 1)).xyz-c;
      fy += textureOffset(map, p, ivec2( 0,-1)).xyz-c;

  return c+D*(fx+fy);
}

void main(){
  fragColor = vec4(diffuse(), 1.0);
}
