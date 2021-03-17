using namespace std;
using namespace glm;

#define PI 3.14159265

double angle(glm::vec2 d){

  if(d.x == 0 && d.y == 0) return 0.0;
  if(d.x == 0 && d.y > 0) return PI/2.0;
  if(d.x == 0 && d.y < 0) return 3.0*PI/2.0;

  double a = 2.0*PI + atan(d.y/d.x);

  if(d.x < 0) a += PI;

  return a;

}

/*
================================================================================
                Individual Segments and Collections (Plates)
================================================================================
*/

struct Litho : Segment {

	Litho(vec2* p):Segment(p){}    //Constructor

  float density = 0.5f;
  float thickness = 0.1f;
  float height = 0.0f;
  float growth = 0.0f;

};

struct Plate {

  Plate(vec2 p):pos{p}{}

  vector<Litho*> seg;

  vec2 pos;
  vec2 speed = vec2(0);
  float rotation = 0.0f;
  float angveloc = 0.0f;
  float mass = 0.0f;
  float inertia = 0.0f;

  //Parameters
  float convection = 150.0f;
  float growth = 0.01f;

  void recenter(){

    pos = vec2(0);
    for(auto&s: seg){
      pos     += *s->pos;
      mass    += s->density*s->thickness;
      inertia += pow(length(pos-*(s->pos)),2)*(s->density*s->thickness);
    }
    pos /= (float)seg.size();

  }
  void update(Cluster<Litho>& clus, double* hm);

};

void Plate::update(Cluster<Litho>& cluster, double* hm){

  //Collide

  for(auto&s: seg){

    int csind = cluster.sample(*(s->pos));

    const int n = 12;
    for(int j = 0; j < n; j++){

      vec2 scan = *(s->pos);
      scan += SIZE*R/12.0f*vec2(cos((float)j/(float)n*2.0f*PI), sin((float)j/(float)n*2.0f*PI));

      if( scan.x >= SIZE || scan.x < 0 ||
          scan.y >= SIZE || scan.y < 0) continue;

      int segind = cluster.sample(scan);

      if(segind < 0) continue;        //Non-Index (Blank Space)
      if(segind == csind) continue;   //Same Segment

      //Add Torque

      s->alive = false;               //Collision
      break;

    }
  }

  //Grow

  for(auto&s: seg){

    if(!s->alive) continue;

    float G = growth*(1.0f-s->thickness);
    ivec2 ip = *(s->pos);

    float nd = 1.0f-hm[ip.y+ip.x*SIZE]; //Hotter = Less Dense

    s->density = s->density*(s->thickness+nd*G/s->density)/(G+s->thickness);
    s->thickness += G;

    s->growth = s->thickness*(1-s->density) - s->height;
    s->height = s->thickness*(1-s->density);

  }

  //Convect

  vec2 acc = vec2(0);
  float torque = 0.0f;

  const function<vec2(ivec2, double*)> force = [](ivec2 i, double* ff){

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

  };

  for(auto&s: seg){

    vec2 f = force(*(s->pos), hm);
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

  //Cleanup

  bool erased = false;
  for(int j = 0; j < seg.size(); j++)
    if(!seg[j]->alive){
      seg.erase(seg.begin()+j);
      j--;
      erased = true;
    }
  if(erased) recenter();

}

/*
================================================================================
                          Main World Container
================================================================================
*/

class World {
public:

  World(int _SEED){

    SEED = _SEED;
    srand(SEED);

    std::cout<<"SEED: "<<SEED<<std::endl;

    perlin.SetOctaveCount(8);
    perlin.SetFrequency(2.0);
    perlin.SetPersistence(0.5);

    heatmap = new double[SIZE*SIZE];
    heightmap = new double[SIZE*SIZE];
    tmpmap = new int[SIZE*SIZE];

    initialize();

    depthmap = new Billboard(SIZE, SIZE);

  }

  ~World(){

    delete[] heatmap;
    delete[] heightmap;
    delete[] tmpmap;

    delete depthmap;
    delete heatA;
    delete heatB;
    delete heightA;
    delete heightB;

  }

  int SEED = 0;
  const glm::vec2 dim = glm::vec2(SIZE, SIZE);
  const float scale = 25.0f;
  noise::module::Perlin perlin;

  vector<Plate> plates;     //Plate Storage
  Cluster<Litho> cluster;   //Cluster Object

  double* heatmap;          //Additional Data
  double* heightmap;
  int* tmpmap;

  Billboard* depthmap;
  Billboard* heatA;
  Billboard* heatB;
  Billboard* heightA;
  Billboard* heightB;

  void initialize();
  void update();

  void subduct(Shader* diffusion, Shader* subduction, Square2D* flat, int n);
  void sediment(Shader* convection, Shader* cascasding, Square2D* flat, int n);

};

/*
================================================================================
                        Initialization and Updating
================================================================================
*/

void World::initialize(){

  //Initialize Heatmap

  float min = 1.0;
  float max = -1.0;
  for(unsigned int i = 0; i < SIZE; i++){
    for(unsigned int j = 0; j < SIZE; j++){
      heatmap[j+i*SIZE] = perlin.GetValue((float)i/(float)SIZE, (float)j/(float)SIZE, SEED);
      if(heatmap[j+i*SIZE] > max) max = heatmap[j+i*SIZE];
      if(heatmap[j+i*SIZE] < min) min = heatmap[j+i*SIZE];
    }
  }

  for(int i = 0; i < SIZE; i++)
    for(int j = 0; j < SIZE; j++)
      heatmap[j+i*SIZE] = (heatmap[j+i*SIZE] - min)/(max-min);

  //Initialize Heightmap

  for(unsigned int i = 0; i < SIZE*SIZE; i++)
    heightmap[i] = 0.31;

  ////Initialize

  //Construct a billboard, using a texture generated from the raw data
  heatA = new Billboard(image::make([&](int i){
    return mix(vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 0.0, 0.0, 1.0), heatmap[i]);
  }, vec2(SIZE, SIZE)));
  heatB = new Billboard(image::make([&](int i){
    return mix(vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 0.0, 0.0, 1.0), heatmap[i]);
  }, vec2(SIZE, SIZE)));

  //Construct a billboard, using a texture generated from the raw data
  heightA = new Billboard(image::make([&](int i){
    return mix(vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0), heightmap[i]);
  }, vec2(SIZE, SIZE)));
  heightB = new Billboard(image::make([&](int i){
    return mix(vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0), heightmap[i]);
  }, vec2(SIZE, SIZE)));

  //Generate Plates
  for(int i = 0; i < nplates; i++)
    plates.emplace_back(Plate(vec2(rand()%SIZE,rand()%SIZE)));

  //Create Relevant clus.segs
  for(auto s: cluster.segs){

    float dist = SIZE*SIZE;
    Plate* nearest;

    for(auto&p: plates){
      if( glm::length(p.pos-*(s->pos)) < dist ){
        dist = glm::length(p.pos-*(s->pos));
        nearest = &p;
      }
    }
    nearest->seg.push_back(s);

  }

  for(auto&p: plates)
    p.recenter();

}

void World::update(){

  // Remove Empty Plates

  for(int i = 0; i < plates.size(); i++){
    if(plates[i].seg.size() == 0){
      plates.erase(plates.begin()+i);
      i--;
      continue;
    }
  }

  // Fill Gaps

  for(auto&p: plates){

    for(auto&s: p.seg){

      float angle = (float)(rand()%100)/100.0f*2.0f*PI;
      vec2 scan = *(s->pos);
      scan += SIZE*R/2.0f*vec2(cos(angle), sin(angle));

      if(scan.x < 0 || scan.x >= SIZE ||
      scan.y < 0 || scan.y >= SIZE) continue;

      //Compute Color at Scan
      //Index of the current guy in general
      int csind = cluster.sample(scan);
      if(csind < 0){

        cluster.points.push_back(scan);
        p.seg.push_back(cluster.add(cluster.points.back()));
        cluster.reassign();
        p.recenter();

        break;

      }
    }
  }

  cluster.update();

}

/*
================================================================================
                            Surface Level Effects
================================================================================
*/

void World::subduct(Shader* diffusion, Shader* subduction, Square2D* flat, int n){

  //Prepare Buffers

  vector<int> colliding;
  for(int i = 0; i < cluster.segs.size(); i++){
    if(!cluster.segs[i]->alive) colliding.push_back(1);
    else colliding.push_back(0);
  }

  //Add SSBO to Shader

  subduction->buffer("colliding", colliding);

  //Execute Alternating Render Pass

  for(int i = 0; i < n; i++){

    heatB->target(false); //No-Clear Target
    diffusion->use();
    diffusion->uniform("D", 0.1f);
    diffusion->uniform("model", mat4(1));
    diffusion->texture("map", heatA->texture);
    flat->render();

    heatA->target(false); //No-Clear Target
    subduction->use();
    subduction->uniform("model", mat4(1));
    subduction->texture("map", heatB->texture);
    subduction->texture("cluster", cluster.target->texture);
    flat->render();

  }

  //Sample heatA into tmpmap

  heatA->sample<int>(tmpmap, vec2(0), dim, GL_COLOR_ATTACHMENT0, GL_RGBA);
  for(int i = 0; i < dim.x*dim.y; i++)
    heatmap[i] = color::i2rgba(tmpmap[i]).r/255.0f;

}

void World::sediment(Shader* convection, Shader* cascading, Square2D* flat, int n){

  //Prepare Buffers

  std::vector<vec2> speed;
  std::vector<float> height;
  std::vector<int> alive;
  for(int i = 0; i < cluster.segs.size(); i++){
    speed.push_back(cluster.segs[i]->speed);
    height.push_back(cluster.segs[i]->height);
    alive.push_back(cluster.segs[i]->alive);
  }

  //Add SSBO to Shader

  cascading->buffer("height", height);
  convection->buffer("speed", speed);
  convection->buffer("alive", alive);

  //Execute Alternating Render Pass

  for(int i = 0; i < n; i++){

    heightB->target(false); //No-Clear Target
    convection->use();
    convection->uniform("model", mat4(1));
    convection->texture("map", heightA->texture);
    convection->texture("cluster", cluster.target->texture);
    flat->render();

    heightA->target(false); //No-Clear Target
    cascading->use();
    cascading->uniform("model", mat4(1));
    cascading->texture("map", heightB->texture);
    cascading->texture("cluster", cluster.target->texture);
    flat->render();

  }

  //Add Sedimentation Offset to Heightmap

  heightA->sample<int>(tmpmap, vec2(0), dim, GL_COLOR_ATTACHMENT0, GL_RGBA);
  vec4 col;
  for(int i = 0; i < dim.x*dim.y; i++){
    col = color::i2rgba(tmpmap[i])/255.0f/255.0f/255.0f;
    heightmap[i] = (col.x+col.y*256+col.z*256*256);
  }

}
