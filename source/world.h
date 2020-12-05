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

#include "plate.h"

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

    perlin.SetOctaveCount(8);
    perlin.SetFrequency(1.0);
    perlin.SetPersistence(0.5);

    initialize();

    clustering = new Billboard(SIZE, SIZE);
    depthmap = new Billboard(SIZE, SIZE);
    heightA = new Billboard(SIZE, SIZE);
    heightB = new Billboard(SIZE, SIZE);

  }

  ~World(){
    delete clustering;
    delete depthmap;
    delete heatA;
    delete heatB;
    delete heightA;
    delete heightB;

    for(int i = 0; i < segments.size(); i++)
      delete segments[i];

    delete[] clustermap;
    delete[] heatmap;
    delete[] heightmap;
  }

  //General Information
  int SEED = 0;
  const glm::vec2 dim = glm::vec2(SIZE, SIZE);
  const float scale = 25.0f;
  noise::module::Perlin perlin;
  float dt = 0.02f;

  double* heatmap;
  double* heightmap;
  int* clustermap;
  int* tmpmap;

  //Plate Centroids
  vector<vec2> centroids;  //Raw Position Buffer
  vector<Litho*> segments; //Segment Pointer Buffer
  vector<Plate> plates;        //Additional Data
  const int nplates = 24;

  Billboard* clustering;
  Billboard* depthmap;
  Billboard* heatA;
  Billboard* heatB;
  Billboard* heightA;
  Billboard* heightB;

  void initialize();
  void drift();
  void cluster(Shader* voronoi, Instance* inst);
  void diffuse(Shader* diffusion, Shader* subduction, Square2D* flat);
  void addRock(Shader* diffusoin, Shader* sedimentation, Square2D* flat);
  void update(Instance* inst);

  void addNode(glm::vec2 pos);
  void delNode(int ind);

  void diffuse(float D, float dt);

};

void World::initialize(){

  heatmap = new double[SIZE*SIZE];
  clustermap = new int[SIZE*SIZE];
  heightmap = new double[SIZE*SIZE];
  tmpmap = new int[SIZE*SIZE];

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

/*
  //Generate Randomized Height
  min = 1.0;
  max = -1.0;
  for(unsigned int i = 0; i < SIZE; i++){
    for(unsigned int j = 0; j < SIZE; j++){
      heightmap[j+i*SIZE] = perlin.GetValue((float)i/(float)SIZE, (float)j/(float)SIZE, SEED);
      if(heightmap[j+i*SIZE] > max) max = heightmap[j+i*SIZE];
      if(heightmap[j+i*SIZE] < min) min = heightmap[j+i*SIZE];
    }
  }

  //Normalize Heightmap
  for(unsigned int i = 0; i < SIZE; i++)
    for(unsigned int j = 0; j < SIZE; j++)
      heightmap[j+i*SIZE] = (heightmap[j+i*SIZE] - min)/(max-min);
*/

 for(unsigned int i = 0; i < SIZE*SIZE; i++)
    heightmap[i] = 0.0;

  //Construct a billboard, using a texture generated from the raw data
  heatA = new Billboard(image::make<double>(vec2(SIZE, SIZE), heatmap, [](double t){
    return mix(vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 0.0, 0.0, 1.0), t);
  }));
  heatB = new Billboard(image::make<double>(vec2(SIZE, SIZE), heatmap, [](double t){
    return mix(vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 0.0, 0.0, 1.0), t);
  }));

  //Construct a billboard, using a texture generated from the raw data
  heightA = new Billboard(image::make<double>(vec2(SIZE, SIZE), heightmap, [](double t){
    return mix(vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0), t);
  }));
  heightB = new Billboard(image::make<double>(vec2(SIZE, SIZE), heightmap, [](double t){
    return mix(vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0), t);
  }));

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

    newseg->plate = nearest;
    nearest->seg.push_back(newseg);

  }

  for(auto&p: plates) p.recenter();

}

/*
================================================================================
                            Plate Dynamics
================================================================================
*/

void World::cluster(Shader* voronoi, Instance* inst){

  clustering->target(glm::vec3(1));
  voronoi->use();
  voronoi->uniform("R", R);
  voronoi->uniform("depthmap", false);
  inst->render();

  clustering->sample<int>(clustermap, vec2(0), dim, GL_COLOR_ATTACHMENT0, GL_RGBA);

}

void World::drift(){

  for(auto& p: plates){
    p.collide(clustermap, centroids, segments);
    p.convect(heatmap, segments);
    p.grow(heatmap);
  }

}

void World::diffuse(Shader* diffusion, Shader* subduction, Square2D* flat){

  std::vector<int> colliding;
  for(int i = 0; i < segments.size(); i++){
    if(segments[i]->colliding) colliding.push_back(1);
    else colliding.push_back(0);
  }
  subduction->buffer("colliding", colliding);

  for(int i = 0; i < 50; i++){

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
    subduction->texture("cluster", clustering->texture);
    flat->render();

  }

  //Sample heatA into tmpmap
  heatA->sample<int>(tmpmap, vec2(0), dim, GL_COLOR_ATTACHMENT0, GL_RGBA);
  for(int i = 0; i < dim.x*dim.y; i++)
    heatmap[i] = color::i2rgba(tmpmap[i]).r/255.0f;

}

bool second = true;
void World::addRock(Shader* convection, Shader* cascading, Square2D* flat){

  std::vector<vec2> speed;
  for(int i = 0; i < segments.size(); i++){
    speed.push_back(segments[i]->speed);
  }
  convection->buffer("speed", speed);

  if(second){
    heightB->target(vec3(0.5)); //Clear the Height Billboard to Black
    heightA->target(vec3(0.5)); //Clear the Height Billboard to Black
    second = false;
  }

  for(int i = 0; i < 5; i++){

    heightB->target(false); //No-Clear Target
    convection->use();
    convection->uniform("model", mat4(1));
    convection->texture("map", heightA->texture);
    convection->texture("cluster", clustering->texture);
    flat->render();

    heightA->target(false); //No-Clear Target
    cascading->use();
    cascading->uniform("model", mat4(1));
    cascading->texture("map", heightB->texture);
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


void World::diffuse(float D, float dt){

};


/*
================================================================================
                            Add / Remove Centroids
================================================================================
*/

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

  newseg->plate = nearest;
	nearest->seg.push_back(newseg);

	for(int i = 0; i < segments.size(); i++)
		segments[i]->pos = &centroids[i]; //Update position as well

}

void World::delNode(int ind){

  centroids.erase(centroids.begin()+ind);
  delete segments[ind];
  segments.erase(segments.begin()+ind);

	for(int i = 0; i < segments.size(); i++)
		segments[i]->pos = &centroids[i]; //Update position as well

}

void World::update(Instance* inst){

  for(int i = 0; i < plates.size(); i++){

    //Remove Empty Plates
    if(plates[i].seg.size() == 0){
      plates.erase(plates.begin()+i);
      i--;
      continue;
    }

    //Remove Colliding Segment References
    bool erased = false;
    for(int j = 0; j < plates[i].seg.size(); j++){
      if(plates[i].seg[j]->colliding){
        plates[i].seg.erase(plates[i].seg.begin()+j);
        j--;
        erased = true;
      }
    }
    if(erased) plates[i].recenter();

  }

  //Remove Colliding Segments
  for(int i = 0; i < segments.size(); i++){
    if(segments[i]->colliding){
      delNode(i);
      i--;
    }
  }

  //This is Inefficient and Requires a better blue noise sampler
  //Add New Nodes
  for(auto&p: plates){
    for(auto&s: p.seg){

      float angle = (float)(rand()%100)/100.0f*2.0f*PI;
      vec2 scan = *(s->pos);
      scan += 256.0f*R/2.0f*vec2(cos(angle), sin(angle));

      if(scan.x < 0 || scan.x >= SIZE ||
      scan.y < 0 || scan.y >= SIZE) continue;

      //Compute Color at Scan
      //Index of the current guy in general
      int cmind = (int)scan.y*SIZE+(int)scan.x;
      vec4 col = color::i2rgba(clustermap[cmind]);

      if(col == vec4(255)){
        addNode(scan);
        p.recenter();
        break;
      }

    }
  }

  inst->updateBuffer(centroids, 0);

}

/*
================================================================================
                                  Rendering
================================================================================
*/

int first = 2;

std::function<void(Model* m, World* w)> tectonicmesh = [](Model* m, World* w){

  if(first) m->indices.clear();
  m->positions.clear();
  m->normals.clear();
  m->colors.clear();

  //Loop over all positions and add the triangles!
  for(int i = 0; i < w->dim.x-1; i++){
    for(int j = 0; j < w->dim.y-1; j++){

      //Get Index
      int ind = i*w->dim.y+j;

      glm::vec4 col = color::i2rgba(w->clustermap[(int)(i*w->dim.y+j)]);
      int aind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(w->clustermap[(int)(i*w->dim.y+j+1)]);
      int bind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(w->clustermap[(int)((i+1)*w->dim.y+j)]);
      int cind = col.x + col.y*256 + col.z*256*256;
      col = color::i2rgba(w->clustermap[(int)((i+1)*w->dim.y+j+1)]);
      int dind = col.x + col.y*256 + col.z*256*256;

      //Add to Position Vector
      glm::vec3 a, b, c, d;

      a = glm::vec3(i  , 0.0, j  );
      b = glm::vec3(i  , 0.0, j+1);
      c = glm::vec3(i+1, 0.0, j  );
      d = glm::vec3(i+1, 0.0, j+1);

      //vec4 stonecolor = glm::vec4(0.8,0.8,0.8,1.0);
      vec4 stonecolor = glm::vec4(0.6,0.68,0.45,1.0);
      vec4 collidecolor = glm::vec4(0.7,0.64,0.52,1.0);
      vec4 magmacolor = glm::vec4(0.84,0.17,0.05,1.0);

      vec4 watercolor = glm::vec4(0.5,0.64,0.87,1.0);
      vec4 earthcolor = glm::vec4(0.68,0.7,0.62,1.0);

      if(viewplates){
        if( aind < w->centroids.size() )
          a += glm::vec3(0, w->scale*w->segments[aind]->height, 0);
        if( bind < w->centroids.size() )
          b += glm::vec3(0, w->scale*w->segments[bind]->height, 0);
        if( cind < w->centroids.size() )
          c += glm::vec3(0, w->scale*w->segments[cind]->height, 0);
        if( dind < w->centroids.size() )
          d += glm::vec3(0, w->scale*w->segments[dind]->height, 0);
      }
      if(!viewplates){

/*
        if( aind < w->centroids.size() )
          a += glm::vec3(0, w->scale*w->segments[aind]->height, 0);
        if( bind < w->centroids.size() )
          b += glm::vec3(0, w->scale*w->segments[bind]->height, 0);
        if( cind < w->centroids.size() )
          c += glm::vec3(0, w->scale*w->segments[cind]->height, 0);
        if( dind < w->centroids.size() )
          d += glm::vec3(0, w->scale*w->segments[dind]->height, 0);
*/

        a += glm::vec3(0, w->scale*w->heightmap[i*(int)w->dim.y+j], 0);
        b += glm::vec3(0, w->scale*w->heightmap[i*(int)w->dim.y+j+1], 0);
        c += glm::vec3(0, w->scale*w->heightmap[(i+1)*(int)w->dim.y+j], 0);
        d += glm::vec3(0, w->scale*w->heightmap[(i+1)*(int)w->dim.y+j+1], 0);

        if(a.y < w->scale*sealevel) a.y = w->scale*sealevel;
        if(b.y < w->scale*sealevel) b.y = w->scale*sealevel;
        if(c.y < w->scale*sealevel) c.y = w->scale*sealevel;
        if(d.y < w->scale*sealevel) d.y = w->scale*sealevel;
      }

      std::function<void(std::vector<GLfloat>&, vec3&)> add3 =
      [](std::vector<GLfloat>& v, vec3& a){
        v.push_back(a.x);
        v.push_back(a.y);
        v.push_back(a.z);
      };

      std::function<void(std::vector<GLfloat>&, vec4&)> add4 =
      [](std::vector<GLfloat>& v, vec4& a){
        v.push_back(a.x);
        v.push_back(a.y);
        v.push_back(a.z);
        v.push_back(a.w);
      };

      //UPPER TRIANGLE

      //Add Indices
      if(first){
        m->indices.push_back(m->positions.size()/3+0);
        m->indices.push_back(m->positions.size()/3+1);
        m->indices.push_back(m->positions.size()/3+2);
      }

      add3(m->positions,a);
      add3(m->positions,b);
      add3(m->positions,c);

      if(viewplates){
        if(aind < w->segments.size()){
          stonecolor = mix(magmacolor, collidecolor, w->segments[aind]->thickness);
          add4(m->colors,stonecolor);
        }
        else add4(m->colors,magmacolor);

        if(bind < w->segments.size()){
          stonecolor = mix(magmacolor, collidecolor, w->segments[bind]->thickness);
          add4(m->colors,stonecolor);
        }
        else add4(m->colors,magmacolor);


        if(cind < w->segments.size()){
          stonecolor = mix(magmacolor, collidecolor, w->segments[cind]->thickness);
          add4(m->colors,stonecolor);
        }
        else add4(m->colors,magmacolor);
      }
      else{
        if(a.y > w->scale*sealevel) add4(m->colors, stonecolor);
        else add4(m->colors, watercolor);
        if(b.y > w->scale*sealevel) add4(m->colors, stonecolor);
        else add4(m->colors, watercolor);
        if(c.y > w->scale*sealevel) add4(m->colors, stonecolor);
        else add4(m->colors, watercolor);
      }

      glm::vec3 n1 = -1.0f*glm::normalize(glm::cross(a-b, c-b));
      for(int i = 0; i < 3; i++)
        add3(m->normals,n1);


      if(first){
        m->indices.push_back(m->positions.size()/3+0);
        m->indices.push_back(m->positions.size()/3+1);
        m->indices.push_back(m->positions.size()/3+2);
      }

      add3(m->positions,d);
      add3(m->positions,c);
      add3(m->positions,b);

      if(viewplates){
        if(dind < w->segments.size()){
          stonecolor = mix(magmacolor, collidecolor, w->segments[dind]->thickness);
          add4(m->colors,stonecolor);
        }
        else add4(m->colors,magmacolor);


        if(cind < w->segments.size()){
          stonecolor = mix(magmacolor, collidecolor, w->segments[cind]->thickness);
          add4(m->colors,stonecolor);
        }
        else add4(m->colors,magmacolor);


        if(bind < w->segments.size()){
          stonecolor = mix(magmacolor, collidecolor, w->segments[bind]->thickness);
          add4(m->colors,stonecolor);
        }
        else add4(m->colors,magmacolor);
      }
      else{
        if(d.y > w->scale*sealevel) add4(m->colors, stonecolor);
        else add4(m->colors, watercolor);
        if(c.y > w->scale*sealevel) add4(m->colors, stonecolor);
        else add4(m->colors, watercolor);
        if(b.y > w->scale*sealevel) add4(m->colors, stonecolor);
        else add4(m->colors, watercolor);
      }

      glm::vec3 n2 = -1.0f*glm::normalize(glm::cross(d-c, b-c));
      for(int i = 0; i < 3; i++)
        add3(m->normals, n2);


    }
  }

  first--;

};
