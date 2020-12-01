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

vec4 cold = vec4(0,0,0,1);
vec4 warm = vec4(1,0,0,1);

void main(){

  //Base Color
  fragColor = texture(map, ex_Tex);

  //Plate Separation
  if(texture(cluster,ex_Tex) == vec4(1.0))
    fragColor = mix(fragColor, cold, 0.2);

  //Plate Collision
  int i = index(texture(cluster, ex_Tex).rgb);
  if(c[i] == 1)
    fragColor = mix(fragColor, warm, 0.2);

}
