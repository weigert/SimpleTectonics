/*
================================================================================
                            Clustered Convection
================================================================================
*/

struct Segment{
public:

  Segment(vec2* p):pos{p}{}

  vec2* pos;
  vec2 speed = vec2(0);
  bool colliding = false;
  int collider = 0;

  float density = 0.5f;
  float thickness = 0.1f;
  float height = 0.0f;
  float growth = 0.0f;

  vec2 force(double* ff){
    const ivec2 i = *pos;
    float fx, fy = 0.0f;

    if(i.x > 0 && i.x < SIZE-1 && i.y > 0 && i.y < SIZE-1){
      fx = -(ff[(i.x+1)*SIZE+i.y] - ff[(i.x-1)*SIZE+i.y])/2.0f;
      fy = -(ff[i.x*SIZE+i.y+1] - ff[i.x*SIZE+i.y-1])/2.0f;
    }

    if(i.x <= 0) fx = 0.1f;
    else if(i.x >= SIZE-1) fx = -0.1f;

    if(i.y <= 0) fy = 0.1;
    else if(i.y >= SIZE-1) fy = -0.1f;

    return vec2(fx,fy);
  }

};

struct Plate {

  Plate(vec2 p):pos{p}{}

  vector<Segment*> seg;

  vec2 pos;
  vec2 speed = vec2(0);
  float rotation = 0.0f;
  float angveloc = 0.0f;
  float mass = 0.0f;
  float inertia = 0.0f;

  //Parameters
  float convection = 150.0f;
  float growth = 0.01f;

  void recenter();
  void convect(double* hm);
  void collide(int* c);
  void grow(double* hm);

};

void Plate::recenter(){

  pos = vec2(0);
  for(auto&s: seg){
    pos     += *s->pos;
    mass    += s->density*s->thickness;
    inertia += pow(length(pos-*(s->pos)),2)*(s->density*s->thickness);
  }
  pos /= (float)seg.size();

}

void Plate::convect(double* hm){

  vec2 acc = vec2(0);
  float torque = 0.0f;

  for(auto&s: seg){

    vec2 f = s->force(hm);
    vec2 dir = *(s->pos)-pos;

    acc += convection*f;
    torque += convection*length(dir)*length(f)*sin(angle(f)-angle(dir));

  }

  speed    += DT*acc/mass;
  angveloc += DT*torque/inertia;
  pos      += DT*speed;
  rotation += DT*angveloc;

  if(rotation > 2*PI) rotation -= 2*PI;
  if(rotation < 0) rotation += 2*PI;

  for(auto&s: seg){

    vec2 dir = *(s->pos) - (pos - DT*speed);
    float _angle = angle(dir) -  (rotation - DT*angveloc);

    s->speed = (pos + length(dir)*vec2(cos(rotation+_angle),sin(rotation+_angle)))-*(s->pos);
    *(s->pos) = pos + length(dir)*vec2(cos(rotation+_angle),sin(rotation+_angle));

  }
}

void Plate::collide(int* c){

  for(auto&s: seg){

    vec2 scan = *(s->pos);

    int cmind = (int)scan.y*SIZE+(int)scan.x;
    vec4 col = color::i2rgba(c[cmind]);
    int csind = (int)col.x + (int)col.y*256 + (int)col.z*256*256;

    const int n = 12;
    for(int j = 0; j < n; j++){

      scan = *(s->pos);
      scan += SIZE*R/12.0f*vec2(cos((float)j/(float)n*2.0f*PI), sin((float)j/(float)n*2.0f*PI));

      if( scan.x >= SIZE || scan.x < 0 ||
          scan.y >= SIZE || scan.y < 0) continue;

      int mapind = (int)scan.y*SIZE+(int)scan.x;
      col = color::i2rgba(c[mapind]);
      int segind = (int)col.x + (int)col.y*256 + (int)col.z*256*256;

      if(col.z == 255) continue;    //Not Colliding
      if(segind == csind) continue; //Not Colliding
      s->colliding = true;          //Colliding
      s->collider = segind;
      break;

    }
  }
}

void Plate::grow(double* hm){

  for(auto&s: seg){

    if(s->colliding) continue;

    float G = growth*(1.0-s->thickness);
    ivec2 ip = *(s->pos);

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
