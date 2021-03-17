/*
================================================================================
                                Rendering Stuff
================================================================================
*/

bool animate = false;
bool viewmap = true;
bool viewplates = true;

float sealevel = 0.32;
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


/*
================================================================================
                                  Rendering
================================================================================
*/

std::function<void(Model* m, World* w)> tectonicmesh = [](Model* m, World* w){

  m->indices.clear();
  m->positions.clear();
  m->normals.clear();
  m->colors.clear();

  //Loop over all positions and add the triangles!
  for(int i = 0; i < w->dim.x-1; i++){
    for(int j = 0; j < w->dim.y-1; j++){

      //Get Index
      int ind = i*w->dim.y+j;

      glm::vec4 col = color::i2rgba(w->cluster.indexmap[(int)(i*w->dim.y+j)]);
      int aind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(w->cluster.indexmap[(int)(i*w->dim.y+j+1)]);
      int bind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(w->cluster.indexmap[(int)((i+1)*w->dim.y+j)]);
      int cind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(w->cluster.indexmap[(int)((i+1)*w->dim.y+j+1)]);
      int dind = col.x + col.y*256 + col.z*256*256;

      //Add to Position Vector
      glm::vec3 a, b, c, d;

      a = glm::vec3(i  , 0.0, j  );
      b = glm::vec3(i  , 0.0, j+1);
      c = glm::vec3(i+1, 0.0, j  );
      d = glm::vec3(i+1, 0.0, j+1);

      float tscale = 50.0f;

      if(viewplates){
        if( aind < w->cluster.points.size() )
          a += glm::vec3(0, w->scale*(w->cluster.segs[aind]->height + w->cluster.segs[aind]->plateheight), 0);
        if( bind < w->cluster.points.size() )
          b += glm::vec3(0, w->scale*(w->cluster.segs[bind]->height + w->cluster.segs[bind]->plateheight), 0);
        if( cind < w->cluster.points.size() )
          c += glm::vec3(0, w->scale*(w->cluster.segs[cind]->height + w->cluster.segs[cind]->plateheight), 0);
        if( dind < w->cluster.points.size() )
          d += glm::vec3(0, w->scale*(w->cluster.segs[dind]->height + w->cluster.segs[dind]->plateheight), 0);
      }
      if(!viewplates){

        a += glm::vec3(0, tscale*(w->heightmap[i*(int)w->dim.y+j]), 0);
        b += glm::vec3(0, tscale*(w->heightmap[i*(int)w->dim.y+j+1]), 0);
        c += glm::vec3(0, tscale*(w->heightmap[(i+1)*(int)w->dim.y+j]), 0);
        d += glm::vec3(0, tscale*(w->heightmap[(i+1)*(int)w->dim.y+j+1]), 0);

        if(a.y < tscale*sealevel) a.y = tscale*sealevel;
        if(b.y < tscale*sealevel) b.y = tscale*sealevel;
        if(c.y < tscale*sealevel) c.y = tscale*sealevel;
        if(d.y < tscale*sealevel) d.y = tscale*sealevel;

      }

      //UPPER TRIANGLE

      //Add Indices
      m->indices.push_back(m->positions.size()/3+0);
      m->indices.push_back(m->positions.size()/3+1);
      m->indices.push_back(m->positions.size()/3+2);

      m->add(m->positions,a);
      m->add(m->positions,b);
      m->add(m->positions,c);
      glm::vec3 n1 = -1.0f*glm::normalize(glm::cross(a-b, c-b));
      for(int i = 0; i < 3; i++)
        m->add(m->normals,n1);

      vec4 tmpcol;

      if(viewplates){
        if(aind < w->cluster.segs.size()){
          tmpcol = mix(magmacolor, collidecolor, w->cluster.segs[aind]->thickness);
          m->add(m->colors,tmpcol);
        }
        else m->add(m->colors,magmacolor);

        if(bind < w->cluster.segs.size()){
          tmpcol = mix(magmacolor, collidecolor, w->cluster.segs[bind]->thickness);
          m->add(m->colors,tmpcol);
        }
        else m->add(m->colors,magmacolor);

        if(cind < w->cluster.segs.size()){
          tmpcol = mix(magmacolor, collidecolor, w->cluster.segs[cind]->thickness);
          m->add(m->colors,tmpcol);
        }
        else m->add(m->colors,magmacolor);
      }
      else{
        if(a.y == tscale*sealevel) m->add(m->colors, watercolor);
        else if(n1.y > steepness) m->add(m->colors, earthcolor);
        else m->add(m->colors, stonecolor);
        if(b.y == tscale*sealevel) m->add(m->colors, watercolor);
        else if(n1.y > steepness) m->add(m->colors, earthcolor);
        else m->add(m->colors, stonecolor);
        if(c.y == tscale*sealevel) m->add(m->colors, watercolor);
        else if(n1.y > steepness) m->add(m->colors, earthcolor);
        else m->add(m->colors, stonecolor);
      }

      m->indices.push_back(m->positions.size()/3+0);
      m->indices.push_back(m->positions.size()/3+1);
      m->indices.push_back(m->positions.size()/3+2);

      m->add(m->positions,d);
      m->add(m->positions,c);
      m->add(m->positions,b);
      glm::vec3 n2 = -1.0f*glm::normalize(glm::cross(d-c, b-c));
      for(int i = 0; i < 3; i++)
        m->add(m->normals, n2);

      if(viewplates){
        if(dind < w->cluster.segs.size()){
          tmpcol = mix(magmacolor, collidecolor, w->cluster.segs[dind]->thickness);
          m->add(m->colors,tmpcol);
        }
        else m->add(m->colors,magmacolor);


        if(cind < w->cluster.segs.size()){
          tmpcol = mix(magmacolor, collidecolor, w->cluster.segs[cind]->thickness);
          m->add(m->colors,tmpcol);
        }
        else m->add(m->colors,magmacolor);


        if(bind < w->cluster.segs.size()){
          tmpcol = mix(magmacolor, collidecolor, w->cluster.segs[bind]->thickness);
          m->add(m->colors,tmpcol);
        }
        else m->add(m->colors,magmacolor);
      }
      else{
        if(d.y == tscale*sealevel) m->add(m->colors, watercolor);
        else if(n2.y > steepness) m->add(m->colors, earthcolor);
        else m->add(m->colors, stonecolor);
        if(c.y == tscale*sealevel) m->add(m->colors, watercolor);
        else if(n2.y > steepness) m->add(m->colors, earthcolor);
        else m->add(m->colors, stonecolor);
        if(b.y == tscale*sealevel) m->add(m->colors, watercolor);
        else if(n2.y > steepness) m->add(m->colors, earthcolor);
        else m->add(m->colors, stonecolor);
      }

    }
  }

};
