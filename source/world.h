using namespace std;
using namespace glm;

#define PI 3.14159265

double angle(glm::vec2 d){

  if(d.x == 0 && d.y == 0) return 0.0;
  if(d.x == 0 && d.y > 0) return PI/2.0;
  if(d.x == 0 && d.y < 0) return 3,0*PI/2.0;

  double a = 2.0*PI + atan(d.y/d.x);

  if(d.x < 0) a += PI;

  return a;

}

/*
==================================================
                Convecting Segments
==================================================
*/

struct Litho{
public:

  Litho(vec2* p):pos{p}{}
  Litho(float d, float t, vec2* p):pos{p},thickness{t},density{d}{}

  glm::vec2* pos;
  float density = 0.5f;
  float thickness = 0.0f;
  float height = 0.0f;
  bool colliding = false;

  vec2 force(double* hm){

    ivec2 i = *pos;
    float fx, fy = 0.0f;

    if(i.x > 0 && i.x < SIZE-1 && i.y > 0 && i.y < SIZE-1){
      fx = (hm[(i.x+1)*SIZE+i.y] - hm[(i.x-1)*SIZE+i.y])/2.0f;
      fy = (hm[i.x*SIZE+i.y+1] - hm[i.x*SIZE+i.y-1])/2.0f;
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

  const float dt = 0.02f;
  const float convection = 150.0f;

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

  void convect(double* hm){

    glm::vec2 acc = glm::vec2(0);
    float torque = 0.0f;

    //Compute Acceleration and Torque
    for(int i = 0 ; i < seg.size(); i++){

      glm::vec2 f = seg[i]->force(hm);
      glm::vec2 dir = *(seg[i]->pos)-pos;

      acc += convection*f;
      torque += convection*length(dir)*length(f)*sin(angle(f)-angle(dir));

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

      *(seg[i]->pos) = pos + length(dir)*vec2(cos(rotation+_angle),sin(rotation+_angle));

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

      else s->colliding = true;

      //Compute Buoyancy
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

      const int n = 8;
      for(int j = 0; j < n; j++){

        //Reset and Shift
        scan = *(s->pos);
        scan += 256.0f*R/10.0f*glm::vec2(cos((float)j/(float)n*2.0f*PI), sin((float)j/(float)n*2.0f*PI));

        //Extract Scan Index
        if( scan.x >= SIZE || scan.x < 0 ||
            scan.y >= SIZE || scan.y < 0) continue;

        //Get Scan Index
        int mapind = (int)scan.y*SIZE+(int)scan.x;
        col = color::i2rgba(c[mapind]);
        int segind = (int)col.x + (int)col.y*256 + (int)col.z*256*256;

        //Plate is not colliding
        if(segind == csind) continue;
      //  if(segs[segind]->colliding) continue;

        //Plate has lower density
      //  if(s->density <= segs[segind]->density) continue;

        s->colliding = true;

        //Plate has higher density and is therefore subducted

        seg.erase(seg.begin()+i);
        i--;
        recenter();

        break;

      }

    }

  }

};

/*
==================================================
              Main World Container
==================================================
*/

class World {
public:

  World(int _SEED){

    SEED = _SEED;
    srand(SEED);

    std::cout<<"SEED: "<<SEED<<std::endl;

    perlin.SetOctaveCount(4);
    perlin.SetFrequency(1.0);
    perlin.SetPersistence(0.5);

    initialize();

    clustering = new Billboard(SIZE, SIZE);
    depthmap = new Billboard(SIZE, SIZE);

  }

  ~World(){
    delete clustering;

    for(int i = 0; i < segments.size(); i++)
      delete segments[i];
  }

  //General Information
  int SEED = 0;
  float t = 0.0;
  const glm::vec2 dim = glm::vec2(SIZE, SIZE);
  noise::module::Perlin perlin;

  double heatmap[SIZE*SIZE]; //Raw Pointer Array (lmao)

  //Plate Centroids
  vector<vec2> centroids;  //Raw Position Buffer
  vector<Litho*> segments; //Segment Pointer Buffer

  vector<Plate> plates;        //Additional Data
  const int nplates = 12;

  Billboard* clustering;
  Billboard* depthmap;
  void initialize();
  void drift(Instance* inst);
  void cluster(Shader* voronoi, Instance* inst);

  void addNode(glm::vec2 pos);
  void delNode(int ind);
};

void World::initialize(){

  //Generate Randomized Heat Map
  float min = 1.0;
  float max = -1.0;
  for(unsigned int i = 0; i < SIZE; i++){
    for(unsigned int j = 0; j < SIZE; j++){
      heatmap[j+i*SIZE] = perlin.GetValue((float)i/(float)SIZE, (float)j/(float)SIZE, SEED);
      if(heatmap[j+i*SIZE] > max) max = heatmap[j+i*SIZE];
      if(heatmap[j+i*SIZE] < min) min = heatmap[j+i*SIZE];
    }
  }

  //Normalize Heatmap
  for(unsigned int i = 0; i < SIZE; i++)
    for(unsigned int j = 0; j < SIZE; j++)
      heatmap[j+i*SIZE] = (heatmap[j+i*SIZE] - min)/(max-min);

  //Generate Plates
  for(int i = 0; i < nplates; i++)
    plates.emplace_back(Plate(vec2(rand()%SIZE,rand()%SIZE)));

  //Generate Plate Centroids
  sample::disc(centroids, K, glm::vec2(0), glm::vec2(256));

  //Create Relevant Segments
  for(auto&c: centroids){

    Litho* newseg = new Litho(0.5f, 0.0, &c); //Properly Scaled Position
    segments.push_back(newseg);

    float dist = SIZE*SIZE;
    Plate* nearest;

    for(auto&p: plates){

      if( glm::length(p.pos-c) < dist ){
        dist = glm::length(p.pos-c);
        nearest = &p;
      }

    }

    nearest->seg.push_back(newseg);

  }

  for(auto&p: plates) p.recenter();

}

void World::cluster(Shader* voronoi, Instance* inst){

  clustering->target(glm::vec3(1));
  voronoi->use();
  voronoi->uniform("R", R);
  voronoi->uniform("depthmap", false);
  inst->render();

/*
  depthmap->target(glm::vec3(1));
  voronoi->use();
  voronoi->uniform("R", R);
  voronoi->uniform("depthmap", true);
  inst->render();
*/
}

void World::addNode(glm::vec2 pos){

  centroids.push_back(pos);
	Litho* newseg = new Litho(0.5f, 0.0f, &centroids.back()); //Properly Scaled Position
	segments.push_back(newseg);

	float dist = SIZE*SIZE;
	Plate* nearest;

	for(auto&p: plates){

		if( glm::length(p.pos-pos) < dist ){
			dist = glm::length(p.pos-pos);
			nearest = &p;
		}

	}

	nearest->seg.push_back(newseg);

	for(int i = 0; i < segments.size(); i++)
		segments[i]->pos = &centroids[i]; //Update position as well

	nearest->recenter();

}

void World::delNode(int ind){

  centroids.erase(centroids.begin()+ind);
  segments.erase(segments.begin()+ind);

	for(int i = 0; i < segments.size(); i++)
		segments[i]->pos = &centroids[i]; //Update position as well

}

void World::drift(Instance* inst){

  int* c = clustering->sample<int>(glm::vec2(0), dim, GL_COLOR_ATTACHMENT0, GL_RGBA);

  for(auto& p: plates){
    p.convect(heatmap);
    p.grow(heatmap);
    p.collide(c, centroids, segments);
  }

  //Remove Colliding Guys?
  for(int i = 0; i < segments.size(); i++){
    if(segments[i]->colliding){
  //    segments[i]->thickness -= 0.1*segments[i]->thickness;
  //    if(segments[i]->thickness < 0.01){
        delNode(i);
        i--;
  //    }
    }
  }

  /*
  Note: At least make it grow at the proper density...
  */

  for(auto&p: plates){
    for(auto&s: p.seg){

      float angle = (float)(rand()%100)/100.0f*2.0f*PI;
      vec2 scan = *(s->pos);
      scan += 256.0f*R/2.0f*vec2(cos(angle), sin(angle));

      //Compute Color at Scan
      //Index of the current guy in general
      int cmind = (int)scan.y*SIZE+(int)scan.x;
      vec4 col = color::i2rgba(c[cmind]);

      if(col == vec4(255)){
        addNode(scan);
        break;
      }

    }
  }

  inst->updateBuffer(centroids, 0);
  delete[] c;

}

/*
==================================================
    Construct a mesh from the world!
==================================================
*/

std::function<void(Model* m, World* w)> tectonicmesh = [](Model* m, World* w){

  m->indices.clear();
  m->positions.clear();
  m->normals.clear();
  m->colors.clear();

  int* inds = w->clustering->sample<int>(glm::vec2(0), w->dim, GL_COLOR_ATTACHMENT0, GL_RGBA);

  //Loop over all positions and add the triangles!
  for(int i = 0; i < w->dim.x-1; i++){
    for(int j = 0; j < w->dim.y-1; j++){

      //Get Index
      int ind = i*w->dim.y+j;

      glm::vec4 col = color::i2rgba(inds[(int)(i*w->dim.y+j)]);
      int aind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(inds[(int)(i*w->dim.y+j+1)]);
      int bind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(inds[(int)((i+1)*w->dim.y+j)]);
      int cind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(inds[(int)((i+1)*w->dim.y+j+1)]);
      int dind = col.x + col.y*256 + col.z*256*256;

      //Add to Position Vector
      glm::vec3 a, b, c, d;

      a = glm::vec3(i  , 0.0, j  );
      b = glm::vec3(i  , 0.0, j+1);
      c = glm::vec3(i+1, 0.0, j  );
      d = glm::vec3(i+1, 0.0, j+1);

      if( aind < w->centroids.size() )
        a = glm::vec3(i  , 10.0*w->segments[aind]->height, j  );
      if( bind < w->centroids.size() )
        b = glm::vec3(i  , 10.0*w->segments[bind]->height, j+1);
      if( cind < w->centroids.size() )
        c = glm::vec3(i+1, 10.0*w->segments[cind]->height, j  );
      if( dind < w->centroids.size() )
        d = glm::vec3(i+1, 10.0*w->segments[dind]->height, j+1);

      glm::vec3 stonecolor = glm::vec3(0.8);
      glm::vec3 collidecolor = glm::vec3(0.7,0.64,0.52);
      glm::vec3 magmacolor = glm::vec3(0.84,0.17,0.05);

      //UPPER TRIANGLE

      //Get the Color of the Ground (Water vs. Flat)

      //Add Indices
      m->indices.push_back(m->positions.size()/3+0);
      m->indices.push_back(m->positions.size()/3+1);
      m->indices.push_back(m->positions.size()/3+2);

      m->positions.push_back(a.x);
      m->positions.push_back(a.y);
      m->positions.push_back(a.z);
      m->positions.push_back(b.x);
      m->positions.push_back(b.y);
      m->positions.push_back(b.z);
      m->positions.push_back(c.x);
      m->positions.push_back(c.y);
      m->positions.push_back(c.z);

      if(aind < w->segments.size()){
        stonecolor = mix(magmacolor, collidecolor, w->segments[aind]->thickness);
        m->colors.push_back(stonecolor.x);
        m->colors.push_back(stonecolor.y);
        m->colors.push_back(stonecolor.z);
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      if(bind < w->segments.size()){
        stonecolor = mix(magmacolor, collidecolor, w->segments[bind]->thickness);
        m->colors.push_back(stonecolor.x);
        m->colors.push_back(stonecolor.y);
        m->colors.push_back(stonecolor.z);
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      if(cind < w->segments.size()){
        stonecolor = mix(magmacolor, collidecolor, w->segments[cind]->thickness);
        m->colors.push_back(stonecolor.x);
        m->colors.push_back(stonecolor.y);
        m->colors.push_back(stonecolor.z);
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      glm::vec3 n1 = glm::normalize(glm::cross(a-b, c-b));

      for(int i = 0; i < 3; i++){
        m->normals.push_back(n1.x);
        m->normals.push_back(n1.y);
        m->normals.push_back(n1.z);
      }

      m->indices.push_back(m->positions.size()/3+0);
      m->indices.push_back(m->positions.size()/3+1);
      m->indices.push_back(m->positions.size()/3+2);

      m->positions.push_back(d.x);
      m->positions.push_back(d.y);
      m->positions.push_back(d.z);
      m->positions.push_back(c.x);
      m->positions.push_back(c.y);
      m->positions.push_back(c.z);
      m->positions.push_back(b.x);
      m->positions.push_back(b.y);
      m->positions.push_back(b.z);

      if(dind < w->segments.size()){
        stonecolor = mix(magmacolor, collidecolor, w->segments[dind]->thickness);
        m->colors.push_back(stonecolor.x);
        m->colors.push_back(stonecolor.y);
        m->colors.push_back(stonecolor.z);
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      if(cind < w->segments.size()){
        stonecolor = mix(magmacolor, collidecolor, w->segments[cind]->thickness);
        m->colors.push_back(stonecolor.x);
        m->colors.push_back(stonecolor.y);
        m->colors.push_back(stonecolor.z);
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      if(bind < w->segments.size()){
        stonecolor = mix(magmacolor, collidecolor, w->segments[bind]->thickness);
        m->colors.push_back(stonecolor.x);
        m->colors.push_back(stonecolor.y);
        m->colors.push_back(stonecolor.z);
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      glm::vec3 n2 = glm::normalize(glm::cross(d-c, b-c));

      for(int i = 0; i < 3; i++){
        m->normals.push_back(n2.x);
        m->normals.push_back(n2.y);
        m->normals.push_back(n2.z);
      }

    }
  }

  delete[] inds;
};

std::function<glm::vec4(double)> heatmap = [](double hm){
  return glm::mix(glm::vec4(0.0,0.0,0.0,1.0), glm::vec4(1.0), hm);
};
