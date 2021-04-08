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

  vec2 speed = vec2(0);       //Velocity of Segment
  bool alive = true;          //Segment Status

  float mass = 0.1f;          //Absolute Accumulated Mass
  float thickness = 0.1f;     //Total Thickness of Plate

  //Area is also known
  float density = 0.5f;       //Density of Mass in Plate (Derived)
  float height = 0.0f;        //Protrusion Height (Buoyant) (Derived)
  float growth = 0.0f;        //Current Growth Rate (Mass / Cycle)

  float plateheight = 0.0f;

};

struct Plate {

  Plate(vec2 p):pos{p}{}

  vector<Litho*> seg;

  vec2 pos;
  vec2 speed = vec2(0);
  float rotation = 0.0f;
  float angveloc = 0.0f;

  float mass = 0.0f;
  float area = 0.0f;

  float inertia = 0.0f;
  float height = 0.0f;

  //Parameters
  float convection = 250.0f;
  float growth = 0.05f;

  void recenter(){

    pos = vec2(0);
    height = 0.0f;
    for(auto&s: seg){
      pos     += *s->pos;

      mass    += s->mass;
      area    += s->area;

      inertia += pow(length(pos-*(s->pos)),2)*(s->mass);
      height  += s->height;
    }
    pos /= (float)seg.size();
    height /= (float)seg.size();

  }
  void update(Cluster<Litho>& clus, double* hm);

};

void Plate::update(Cluster<Litho>& cluster, double* hm){

  vec2 acc = vec2(0);
  float torque = 0.0f;

  //Collide

  for(auto&s: seg){

    ivec2 ipos = *(s->pos);

    if( ipos.x >= SIZE || ipos.x < 0 ||
        ipos.y >= SIZE || ipos.y < 0){
          *(s->pos);
          s->alive = false;
          continue;
        }

    int csind = cluster.sample(ipos);

    const int n = 12;
    for(int j = 0; j < n; j++){

      vec2 scan = *(s->pos);
      scan += SIZE*CR*vec2(cos((float)j/(float)n*2.0f*PI), sin((float)j/(float)n*2.0f*PI));

      if( scan.x >= SIZE || scan.x < 0 ||
          scan.y >= SIZE || scan.y < 0) continue;

      int segind = cluster.sample(scan);

      if(segind < 0) continue;        //Non-Index (Blank Space)
      if(segind == csind) continue;   //Same Segment


      //Two Segments are Colliding, Subduce the Denser One
      if(cluster.segs[csind]->height > cluster.segs[segind]->height){

      //  cluster.segs[segind]->thickness += 0.2*cluster.segs[csind]->height;  //Move Mass
        s->alive = false;
        break;

      }

    }
  }

  //Grow

  const std::function<float(float, float)> equDensity = [&](float k, float temp){
    return k*temp/(1.0f+k*temp);
  };

  for(auto&s: seg){

    if(!s->alive) continue;

    ivec2 ip = *(s->pos);
    float nd = hm[ip.y+ip.x*SIZE];              //Heat Value!


    //LINEAR GROWTH RATE [m / s]

    float G = growth*(1.1-nd)*(nd-s->height);
    if(G < 0.0) G *= 0.1;

    //COMPUTE EQUILIBRIUM DENSITY (PER-VOLUME)
    float D = equDensity(1.0f, nd);

    s->mass += s->area*G*D; //m^2 * m / s * kg / m^3 = kg
    s->thickness = s->thickness + G; //New Thickness
    s->density = s->mass/(s->area*s->thickness);

    //Height is Computed
    s->height = s->thickness*(1.0f-s->density);

  }

  //Convect

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

    vec2 f = -1.0f*force(*(s->pos), hm);
    vec2 dir = *(s->pos)-pos;

    acc -= convection*f;
    torque -= convection*length(dir)*length(f)*sin(angle(f)-angle(dir));

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
  const float scale = 50.0f;
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
  void sediment(Shader* cascasding, Square2D* flat, int n);

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

  for(auto&p: plates){

/*
    for(auto&s: p.seg){
      ivec2 ip = *(s->pos);
      if(ip.x < -SIZE || ip.x > 2*SIZE-1 ||
      ip.y < -SIZE || ip.y > 2*SIZE-1) s->alive = false;
    }
    */

    bool erased = false;
    for(int j = 0; j < p.seg.size(); j++)
      if(!p.seg[j]->alive){
        p.seg.erase(p.seg.begin()+j);
        j--;
        erased = true;
      }
    if(erased) p.recenter();

  }

  for(int i = 0; i < plates.size(); i++){
    if(plates[i].seg.size() == 0){
      plates.erase(plates.begin()+i);
      i--;
      continue;
    }
  }

  cluster.remove([](Litho* s){
    return !s->alive;
  });

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
    diffusion->uniform("D", 0.2f);
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

void World::sediment(Shader* cascading, Square2D* flat, int n){

  //Prepare Buffers

  std::vector<float> height;
  for(int i = 0; i < cluster.segs.size(); i++){
    height.push_back(cluster.segs[i]->height+cluster.segs[i]->plateheight);
  }

  //Add SSBO to Shader

  cascading->buffer("height", height);

  //Execute Alternating Render Pass

  heightA->target(false); //No-Clear Target
  cascading->use();
  cascading->uniform("model", mat4(1));
  cascading->uniform("init", true);
  cascading->texture("map", heightB->texture);
  cascading->texture("cluster", cluster.target->texture);
  flat->render();

  for(int i = 0; i < n; i++){

    heightB->target(false); //No-Clear Target
    cascading->use();
    cascading->uniform("model", mat4(1));
    cascading->uniform("init", false);
    cascading->texture("map", heightA->texture);
    cascading->texture("cluster", cluster.target->texture);
    flat->render();

    heightA->target(false); //No-Clear Target
    cascading->use();
    cascading->uniform("model", mat4(1));
    cascading->uniform("init", false);
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
