

/*
================================================================================
                              Convecting Segments
================================================================================
*/

struct Plate;

struct Litho{
public:

  Litho(vec2* p):pos{p}{}
  Litho(float d, float t, vec2* p):pos{p},thickness{t},density{d}{}
  Litho(float d, float t, vec2* p, Plate* pl):pos{p},thickness{t},density{d},plate{pl}{}

  glm::vec2* pos;
  vec2 speed = vec2(0);
  float density = 0.5f;
  float thickness = 0.0f;
  float height = 0.0f;
  bool colliding = false;
  int collider = 0;
  float growth = 0.0f;

  //Parent
  Plate* plate;

  vec2 force(double* hm){

    const ivec2 i = *pos;
    float fx, fy = 0.0f;

    if(i.x > 0 && i.x < SIZE-1 && i.y > 0 && i.y < SIZE-1){
      fx = -(hm[(i.x+1)*SIZE+i.y] - hm[(i.x-1)*SIZE+i.y])/2.0f;
      fy = -(hm[i.x*SIZE+i.y+1] - hm[i.x*SIZE+i.y-1])/2.0f;
    }

    if(i.x <= 0) fx = 0.1f;
    else if(i.x >= SIZE-1) fx = -0.1f;

    if(i.y <= 0) fy = 0.1;
    else if(i.y >= SIZE-1) fy = -0.1f;

    return vec2(fx,fy);

  }

};

/*
==================================================
          Clustered Convection Structure
==================================================
*/


struct Plate {

  Plate(vec2 p):pos{p}{}

  std::vector<Litho*> seg;

  vec2 pos;
  vec2 speed = vec2(0);
  float rotation = 0.0f;
  float angveloc = 0.0f;
  float mass = 0.0f;
  float inertia = 0.0f;

  float dt = 0.02f;
  float convection = 150.0f;

  void recenter(){

      pos = vec2(0);

      for(int i = 0; i < seg.size(); i++){
        pos += vec2((*seg[i]).pos->x, (*seg[i]).pos->y);
      }
      pos /= vec2(seg.size());

      //Mass and Inertia
      mass = seg.size();
      for(int i = 0; i < seg.size(); i++){
        inertia += pow(length(pos-*(seg[i]->pos)),2);
      }

  }

  void convect(double* hm, vector<Litho*>& segments){

    glm::vec2 acc = glm::vec2(0);
    float torque = 0.0f;

    //Compute Acceleration and Torque
    for(int i = 0 ; i < seg.size(); i++){

      glm::vec2 f = seg[i]->force(hm);

      //Collision Force
      if(seg[i]->colliding){

        //Collision Force Direction
        vec2 cf = *seg[i]->pos - *segments[seg[i]->collider]->pos;
        if(length(*seg[i]->pos - pos) > length(*segments[seg[i]->collider]->pos - pos))
          cf *= -1.0f;

        if(length(cf) != 0)
          f += 2.5f*normalize(cf);

      }

      glm::vec2 dir = *(seg[i]->pos)-pos;

      acc += convection*f;
      torque += 5*convection*length(dir)*length(f)*sin(angle(f)-angle(dir));

    }

    //Terminal Acceleration
    acc*= (5.0-length(speed));

    //Change Speed and Angular Velocity
    speed += dt*acc/mass;
    angveloc += dt*torque/inertia;

    //Move Plate
    pos += dt*speed;
    rotation += dt*angveloc;
    if(rotation > 2*PI) rotation -= 2*PI;
    if(rotation < 0) rotation += 2*PI;

    glm::vec2 dir;
    float _angle;

    //Move Segments
    for(int i = 0 ; i < seg.size(); i++){

      dir = *(seg[i]->pos) - (pos - dt*speed);
      _angle = angle(dir) -  (rotation - dt*angveloc);

      seg[i]->speed = length(dir)*vec2(cos(rotation+_angle),sin(rotation+_angle));
      *(seg[i]->pos) = pos + seg[i]->speed;

    }

  }

  void grow(double* hm){

    for(int i = 0 ; i < seg.size(); i++){
      Litho* s = seg[i];

      if(s->colliding) continue;

      glm::ivec2 ip = *(s->pos);

      float G = 0.01f*(1.0-s->thickness);

      if(ip.x >= 0 && ip.x < SIZE-1 &&
      ip.y >= 0 && ip.y < SIZE-1){

        float nd = 1.0-hm[ip.y+ip.x*SIZE]; //Hotter = Less Dense

        s->density = s->density*(s->thickness+nd*G/s->density)/(G+s->thickness);
        s->thickness += G;

      }
      else if(ip.x < -SIZE || ip.x > 2*SIZE-1 ||
      ip.y < -SIZE || ip.y > 2*SIZE-1){
        s->colliding = true;
      }

      //Compute Buoyancy
      s->growth = s->thickness*(1-s->density) - s->height;
      s->height = s->thickness*(1-s->density);

    }
  }

  void collide(int* c, vector<vec2>& centroids, vector<Litho*>& segs){

    for(int i = 0 ; i < seg.size(); i++){

      Litho* s = seg[i];
      s->colliding = false;

      glm::vec2 scan = *(s->pos);

      //Index of the current guy in general
      int cmind = (int)scan.y*SIZE+(int)scan.x;
      glm::vec4 col = color::i2rgba(c[cmind]);
      int csind = (int)col.x + (int)col.y*256 + (int)col.z*256*256;

      const int n = 12;
      for(int j = 0; j < n; j++){

        //Reset and Shift
        scan = *(s->pos);
        scan += 256.0f*R/12.0f*glm::vec2(cos((float)j/(float)n*2.0f*PI), sin((float)j/(float)n*2.0f*PI));

        //Extract Scan Index
        if( scan.x >= SIZE || scan.x < 0 ||
            scan.y >= SIZE || scan.y < 0) continue;

        //Get Scan Index
        int mapind = (int)scan.y*SIZE+(int)scan.x;
        col = color::i2rgba(c[mapind]);

        if(col.z == 255) continue;

        int segind = (int)col.x + (int)col.y*256 + (int)col.z*256*256;

        //Plate is not colliding
        if(segind == csind) continue;

      //  if(s->density < segs[segind]->density) continue;

        s->colliding = true;
        s->collider = segind;

        break;

      }

    }

  }

};
