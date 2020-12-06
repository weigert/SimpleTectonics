#version 330 core
in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;

uniform float D;

vec3 diffuse(float rate){
  vec2 p = ex_Tex;
  vec3 c = textureOffset(map, p, ivec2(0,0)).xyz;

  vec3 diff = textureOffset(map, p, ivec2( 1, 0)).xyz-c;
      diff += textureOffset(map, p, ivec2(-1, 0)).xyz-c;
      diff += textureOffset(map, p, ivec2( 0, 1)).xyz-c;
      diff += textureOffset(map, p, ivec2( 0,-1)).xyz-c;

  return c+rate*diff;
}

void main(){
  fragColor = vec4(diffuse(D), 1.0);
}
