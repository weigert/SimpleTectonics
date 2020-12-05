#version 430 core
in vec2 ex_Tex;
out vec4 fragColor;

uniform sampler2D map;
uniform sampler2D cluster;

layout (std430, binding = 0) buffer colliding {
  int c[];
};

int index(vec3 a){
  a = a*255.0f;
  return int((a.x + a.y*256 + a.z*256*256));
}

float num(vec3 a){
  a = a*255.0f;
  return (a.x + a.y*256 + a.z*256*256);
}

vec3 col(int i){
  float r = ((i >>  0) & 0xff)/255.0f;
  float g = ((i >>  8) & 0xff)/255.0f;
  float b = ((i >> 16) & 0xff)/255.0f;
  return vec3(r,g,b);
}

vec4 cold = vec4(0,0,0,1);
vec4 warm = vec4(1,1,1,1);

void main(){

  //Base Color
  fragColor = texture(map, ex_Tex);

/*
  //Plate Separation
  if(texture(cluster,ex_Tex) == vec4(1.0))
    fragColor -= vec4(col(100000),0.0);

  //Plate Collision
  int i = index(texture(cluster, ex_Tex).rgb);
  if(c[i] == 1)
    fragColor += vec4(col(100), 0.0);
*/
}
