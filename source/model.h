//GPU Accelerated Voronoi Controller Stuff
bool animate = false;
float K = 512;
float R = 2.0f*sqrt(4.0f/3.14159265f/K);

#define SIZE 256

/*
Rendering Stuff
*/

const int WIDTH = 800;
const int HEIGHT = 800;

float zoom = 0.2;
float zoomInc = 0.005;

//Rotation and View
float rotation = 0.0f;
glm::vec3 cameraPos = glm::vec3(50, 50, 50);
glm::vec3 lookPos = glm::vec3(0, 0, 0);
glm::mat4 camera = glm::lookAt(cameraPos, lookPos, glm::vec3(0,1,0));
glm::mat4 projection = glm::ortho(-(float)WIDTH*zoom, (float)WIDTH*zoom, -(float)HEIGHT*zoom, (float)HEIGHT*zoom, -800.0f, 500.0f);

glm::vec3 viewPos = glm::vec3(SIZE/2.0, 40.0, SIZE/2.0);

//Shader Stuff
float steepness = 0.8;
glm::vec3 flatColor = glm::vec3(0.40, 0.60, 0.25);
glm::vec3 waterColor = glm::vec3(0.17, 0.20, 0.4);
glm::vec3 steepColor = glm::vec3(0.7);
//glm::vec3 steepColor = glm::vec3(0.78, 0.6, 0.168);
//glm::vec3 flatColor = glm::vec3(0.84, 0.65, 0.36);

//Lighting and Shading
glm::vec3 skyCol = glm::vec4(0.64, 0.75, 0.9, 1.0f);
glm::vec3 lightPos = glm::vec3(-100.0f, 100.0f, -150.0f);
glm::vec3 lightCol = glm::vec3(1.0f, 1.0f, 0.9f);
float lightStrength = 1.4;
glm::mat4 depthModelMatrix = glm::mat4(1.0);
glm::mat4 depthProjection = glm::ortho<float>(-300, 300, -300, 300, 0, 800);
glm::mat4 depthCamera = glm::lookAt(lightPos, glm::vec3(0), glm::vec3(0,1,0));
bool viewmap = true;

glm::mat4 biasMatrix = glm::mat4(
    0.5, 0.0, 0.0, 0.0,
    0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.5, 0.0,
    0.5, 0.5, 0.5, 1.0
);


std::function<void()> eventHandler = [&](){

  bool changed = false;

    if(Tiny::event.scroll.posy && zoom <= 0.3){
      zoom /= 0.975;
      projection = glm::ortho(-(float)WIDTH*zoom, (float)WIDTH*zoom, -(float)HEIGHT*zoom, (float)HEIGHT*zoom, -800.0f, 500.0f);
      changed = true;
    }
    else if(Tiny::event.scroll.negy && zoom > 0.005){
      zoom *= 0.975;
      projection = glm::ortho(-(float)WIDTH*zoom, (float)WIDTH*zoom, -(float)HEIGHT*zoom, (float)HEIGHT*zoom, -800.0f, 500.0f);
      changed = true;
    }
    else if(Tiny::event.scroll.negx){
      rotation += 1.5f;
      camera = glm::rotate(camera, glm::radians(1.5f), glm::vec3(0.0f, 1.0f, 0.0f));
      changed = true;
    }
    else if(Tiny::event.scroll.posx){
      rotation -= 1.5f;
      camera = glm::rotate(camera, glm::radians(-1.5f), glm::vec3(0.0f, 1.0f, 0.0f));
      changed = true;
    }

  //Adjust Stuff
  if(changed){
    if(rotation < 0.0) rotation = 360.0 + rotation;
    else if(rotation > 360.0) rotation = rotation - 360.0;
    camera = glm::rotate(glm::lookAt(cameraPos, lookPos, glm::vec3(0,1,0)), glm::radians(rotation), glm::vec3(0,1,0));
  }

  if(!Tiny::event.press.empty()){

    if(Tiny::event.press.back() == SDLK_p)
      animate = !animate;

    if(Tiny::event.press.back() == SDLK_ESCAPE)
      viewmap = !viewmap;

    if(Tiny::event.press.back() == SDLK_SPACE)
      viewPos += glm::vec3(0.0, 1.0, 0.0);

    if(Tiny::event.press.back() == SDLK_c)
      viewPos -= glm::vec3(0.0, 1.0, 0.0);

    if(Tiny::event.press.back() == SDLK_w)
      viewPos -= glm::vec3(1.0, 0.0, 0.0);

    if(Tiny::event.press.back() == SDLK_a)
      viewPos += glm::vec3(0.0, 0.0, 1.0);

    if(Tiny::event.press.back() == SDLK_s)
      viewPos += glm::vec3(1.0, 0.0, 0.0);

    if(Tiny::event.press.back() == SDLK_d)
      viewPos -= glm::vec3(0.0, 0.0, 1.0);

    if(Tiny::event.press.back() == SDLK_UP){
      cameraPos += glm::vec3(0, 5, 0);
      camera = glm::rotate(glm::lookAt(cameraPos, lookPos, glm::vec3(0,1,0)), glm::radians(rotation), glm::vec3(0,1,0));
    }

    if(Tiny::event.press.back() == SDLK_DOWN){
      cameraPos -= glm::vec3(0, 5, 0);
      camera = glm::rotate(glm::lookAt(cameraPos, lookPos, glm::vec3(0,1,0)), glm::radians(rotation), glm::vec3(0,1,0));
    }

  }

};
