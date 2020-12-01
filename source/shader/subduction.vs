#version 330 core
in vec2 in_Quad;
in vec2 in_Tex;
out vec2 ex_Tex;

//Position the Billboard in space!
const mat4 model = mat4(1);

void main(){
  ex_Tex = in_Tex;
  gl_Position = model*vec4(in_Quad, 1.0, 1.0);
}
