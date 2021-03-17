/*
================================================================================
                                Rendering Stuff
================================================================================
*/

bool animate = false;
bool viewmap = true;
bool viewplates = true;

float sealevel = 0.0;
float steepness = 0.95;

float zoom = 0.2;
float zoomInc = 0.005;
float rotation = 180.0f;

glm::vec3 cameraPos = glm::vec3(-50, 50, -50);
glm::vec3 lookPos = glm::vec3(0, 0, 0);
glm::mat4 camera = glm::rotate(glm::lookAt(cameraPos, lookPos, glm::vec3(0,1,0)), glm::radians(rotation), glm::vec3(0.0f, 1.0f, 0.0f));
glm::vec3 viewPos = glm::vec3(SIZE/2.0, 40.0, SIZE/2.0);

glm::mat4 projection = glm::ortho(-(float)WIDTH*zoom, (float)WIDTH*zoom, -(float)HEIGHT*zoom, (float)HEIGHT*zoom, -800.0f, 500.0f);

glm::vec3 skyCol = glm::vec4(0.17, 0.11, 0.18, 1.0f);
glm::vec3 skyBlue = glm::vec4(0.16, 0.14, 0.14, 1.0f);
glm::vec4 collidecolor = glm::vec4(0.7,0.64,0.52,1.0);
glm::vec4 magmacolor = glm::vec4(0.84,0.17,0.05,1.0);
glm::vec4 watercolor = glm::vec4(0.5,0.64,0.87,1.0);
glm::vec4 earthcolor = glm::vec4(0.89,0.78,0.73,1.0);
glm::vec4 stonecolor = glm::vec4(0.77,0.72,0.71,1.0);

glm::vec3 lightPos = glm::vec3(-100.0f, 120.0f, -150.0f);
glm::vec3 lightCol = glm::vec3(1.0f, 1.0f, 0.9f);
float lightStrength = 1.4;

float sb[3] = {skyBlue.x, skyBlue.y, skyBlue.z};
float sc[3] = {stonecolor.x, stonecolor.y, stonecolor.z};
float ec[3] = {earthcolor.x, earthcolor.y, earthcolor.z};

glm::mat4 depthModelMatrix = glm::mat4(1.0);
glm::mat4 depthProjection = glm::ortho<float>(-300, 300, -300, 300, 0, 800);
glm::mat4 depthCamera = glm::lookAt(lightPos, glm::vec3(0), glm::vec3(0,1,0));
glm::mat4 biasMatrix = glm::mat4(
    0.5, 0.0, 0.0, 0.0,
    0.0, 0.5, 0.0, 0.0,
    0.0, 0.0, 0.5, 0.0,
    0.5, 0.5, 0.5, 1.0
);

Handle interfaceFunc = [](){

  ImGui::Text("Simulation Controller");
  ImGui::DragFloat("Sealevel", &sealevel, 0.01, 0, 1);
  ImGui::DragFloat("Steepness", &steepness, 0.01, 0, 1);
  ImGui::ColorEdit3("Sky Color", sb);
  ImGui::ColorEdit3("Ground Color", ec);
  ImGui::ColorEdit3("Stone Color", sc);
  if(ImGui::Button("Set")){
    skyBlue.x = sb[0];
    skyBlue.y = sb[1];
    skyBlue.z = sb[2];
    earthcolor.x = ec[0];
    earthcolor.y = ec[1];
    earthcolor.z = ec[2];
    stonecolor.x = sc[0];
    stonecolor.y = sc[1];
    stonecolor.z = sc[2];
  }

};

std::function<void()> eventHandler = [](){

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

  if(Tiny::event.active[SDLK_SPACE])
    viewPos += glm::vec3(0.0, 1.0, 0.0);
  if(Tiny::event.active[SDLK_c])
    viewPos -= glm::vec3(0.0, 1.0, 0.0);
  if(Tiny::event.active[SDLK_UP]){
    cameraPos += glm::vec3(0, 1, 0);
    camera = glm::rotate(glm::lookAt(cameraPos, lookPos, glm::vec3(0,1,0)), glm::radians(rotation), glm::vec3(0,1,0));
  }
  if(Tiny::event.active[SDLK_DOWN]){
    cameraPos -= glm::vec3(0, 1, 0);
    camera = glm::rotate(glm::lookAt(cameraPos, lookPos, glm::vec3(0,1,0)), glm::radians(rotation), glm::vec3(0,1,0));
  }

  if(Tiny::event.active[SDLK_w])
    viewPos -= glm::vec3(1.0, 0.0, 0.0);
  if(Tiny::event.active[SDLK_a])
    viewPos += glm::vec3(0.0, 0.0, 1.0);
  if(Tiny::event.active[SDLK_s])
    viewPos += glm::vec3(1.0, 0.0, 0.0);
  if(Tiny::event.active[SDLK_d])
    viewPos -= glm::vec3(0.0, 0.0, 1.0);

  if(!Tiny::event.press.empty()){

    if(Tiny::event.press.back() == SDLK_p)
      animate = !animate;

    if(Tiny::event.press.back() == SDLK_m)
      viewmap = !viewmap;

    if(Tiny::event.press.back() == SDLK_RETURN)
      viewplates = !viewplates;


  }

};
