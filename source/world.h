/*
    Plates track their Centroids

    //Plate Movement
    Based on a heat map, we can compute centroid convection
    Convection is proportional to the gradient obvisously.
    The heat map is initialized with noise.
    The centroid convections are used to compute the force and torque of the plate
    we can compute the plate center, direction and rotation
    we can then move the centroids

    //Centroid Height
    Centroids have a thickness and a density.
    Age of a segment dictates its thickness, which grows at a certain rate,
    while the heat dictates an equilibrium density which is the density the plate grows at.
    The centroid thereby stores a certain amount of mass over its area.

    Additionally, the centroid's 3D position is determined by its density and volume.
    It is basically floating, and there is a certain amount above and below the surface.
    Density is initially 0.5, so that it floats half above.

    //Plate Collision
    When plates collide, the less dense plate is pushed below and the more dense plate is pushed above.
    This creates sediment and heat. How much sediment and how much heat? Basically we take the combined
    height and density, compute the overshoot on the top and bottom to compute sedimentation and heating rate.

    When plates separate, new plate is created which has the density given by the heat below.
    That is all.

    //Heat Generation
    Heat can initially be seeded with random noise and we can see how plates move.
    Later heat can be generated randomly or can be affected by the centroid growth.
    When centroids collide, they generate heat by destruction of mass.

    The area of the less dense plate is destroyed proportionally to its submerged density

    //Sediment Transport and Heat Generation
    Computing where and how plates move relative to each other,
    we can compute where to generate heat.

    NOTE:
    As a plate with a high density, i.e. sitting at a ratio other than 1:1 in the ground,
    then as the thickness decreases due to collision, the plate will sink further down faster.
    The thickness of the plate needs to decrease as it subduces.

*/

using namespace std;
using namespace glm;

#define PI 3.14159265

//Plate Centroid
struct Litho{
public:

  Litho(glm::vec2* p){
    pos = p;
  }

  Litho(float d, float t, glm::vec2* p){
    density = d; thickness = t; pos = p;
  }

  glm::vec2* pos;
  float density = 0.5f;
  float thickness = 0.0f;
  float height = 0.0f;
  bool colliding = false;

  glm::vec2 force(double* hm){

    glm::ivec2 i = *pos;
    float fx, fy = 0.0f;

    //Fully in Bound
    if(i.x > 0 && i.x < SIZE-1 && i.y > 0 && i.y < SIZE-1){
      fx = (hm[(i.x+1)*SIZE+i.y] - hm[(i.x-1)*SIZE+i.y])/2.0f;
      fy = (hm[i.x*SIZE+i.y+1] - hm[i.x*SIZE+i.y-1])/2.0f;
    }

    if(i.x <= 0) fx = 0.1f;
    else if(i.x >= SIZE-1) fx = -0.1f;

    if(i.y <= 0) fy = 0.1;
    else if(i.y >= SIZE-1) fy = -0.1f;

    return glm::vec2(fx,fy);

  }

};

double angle(glm::vec2 d){

  if(d.x == 0 && d.y == 0) return 0.0;
  if(d.x == 0 && d.y > 0) return PI/2.0;
  if(d.x == 0 && d.y < 0) return 3,0*PI/2.0;

  double a = 2.0*PI + atan(d.y/d.x);

  if(d.x < 0) a += PI;

  return a;

}

struct Plate {

  std::vector<int> segments;

  glm::vec2 pos;
  glm::vec2 speed = glm::vec2(0);
  float rotation = 0.0f;
  float angveloc = 0.0f;
  float mass = 0.0f;
  float inertia = 0.0f;

  const float dt = 0.02f;
  const float convection = 200.0f;

  void recenter(vector<Litho>& s){
    pos = glm::vec2(0);
    for(auto&i: segments)
      pos += *s[i].pos;
    pos /= glm::vec2(segments.size());

    //Mass and Inertia
    mass = segments.size();
    for(auto&i: segments)
      inertia += pow(length(pos-*s[i].pos),2);
  }

  void convect(double* hm, vector<Litho>& s){

    glm::vec2 acc = glm::vec2(0);
    float torque = 0.0f;

    //Compute Acceleration and Torque
    for(auto&i: segments){

      glm::vec2 f = s[i].force(hm);
      glm::vec2 dir = *s[i].pos-pos;

      acc += convection*f;
      torque += convection*length(dir)*length(f)*sin(angle(f)-angle(dir));

    }

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
    for(auto&i: segments){

      dir = *s[i].pos - (pos - dt*speed);
      _angle = angle(dir) -  (rotation - dt*angveloc);

      *s[i].pos = pos + length(dir)*vec2(cos(rotation+_angle),sin(rotation+_angle));

    }

  }

  void grow(double* hm, std::vector<Litho>& s){

    for(auto&i: segments){

      glm::ivec2 ip = *s[i].pos;

      float G = 0.01f*(1.0-s[i].thickness);

      if(ip.x >= 0 && ip.x < SIZE-1 &&
      ip.y >= 0 && ip.y < SIZE-1){

        float nd = 1.0-hm[ip.y+ip.x*SIZE]; //Hotter = Less Dense

        s[i].density = s[i].density*(s[i].thickness+nd*G/s[i].density)/(G+s[i].thickness);
        s[i].thickness += G;

      }

      //Compute Buoyancy
      s[i].height = s[i].thickness*(1-s[i].density);

    }
  }

  void collide(int* c, std::vector<Litho>& s){

    /*
        Iterate over the segments in your guy, identify ones which have
        an alpha value of over two meaning they overlap with the area of another guy

        Inside a collision radius R/2, we check for other segments by scanning the image.
        Then we can use that to collide.

        Inside a segment, you can't collide with yourself! No need to track.
    */


    for(auto& i: segments){

      glm::vec2 scan;

      s[i].colliding = false;

      int n = 8;
      for(int j = 0; j < n; j++){

        scan = *s[i].pos;
        scan += 256.0f*R/10.0f*glm::vec2(cos((float)j/(float)n*2.0f*PI), sin((float)j/(float)n*2.0f*PI));

        if( scan.x >= SIZE || scan.x < 0 ||
            scan.y >= SIZE || scan.y < 0) continue;

        //std::cout<<i<<": "<<s[i].pos->x<<" "<<s[i].pos->y<<std::endl;

        //Scan Position on Map
        int mapind = (int)scan.y*SIZE+(int)scan.x; //This is flipped for some reason in the texture

        glm::vec4 col = color::i2rgba(c[mapind]);
        int segind = col.x + col.y*256 + col.z*256*256;

        //Collision between segment i and sind (for visualization)
        if(segind == i) continue;

        s[i].colliding = true;

        if(s[i].density > s[segind].density){
          s[i].thickness -= 0.05*s[i].thickness;

          //Now: Remove elements from the vector!

          if(s[i].thickness < 0.0001);
          //Delete the centroid and simultaneously
            //Remove guy
        }
        else{
          s[segind].thickness -= 0.05*s[segind].thickness;
        }

      //  if(s[i].thickness < )



      }




      //std::cout<<sind<<std::endl;

    }

  }

};

//Prepare Noise for jiggling the centroids

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
  }

  //General Information
  int SEED = 0;
  float t = 0.0;
  const glm::vec2 dim = glm::vec2(SIZE, SIZE);
  noise::module::Perlin perlin;

  double heatmap[SIZE*SIZE]; //Raw Pointer Array (lmao)

  //Plate Centroids
  std::vector<glm::vec2> centroids; //Raw Position Buffer
  std::vector<Litho> segments; //Raw Position Buffer

  std::vector<Plate> plates;        //Additional Data
  const int nplates = 12;

  Billboard* clustering;
  Billboard* depthmap;
  void initialize();
  void drift(Instance* inst);
  void cluster(Shader* voronoi, Instance* inst);

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

  for(unsigned int i = 0; i < SIZE; i++){
    for(unsigned int j = 0; j < SIZE; j++){
      heatmap[j+i*SIZE] = (heatmap[j+i*SIZE] - min)/(max-min);
    }
  }

  //Generate Plates
  for(int i = 0; i < nplates; i++){

    Plate plate;
    plate.pos = glm::vec2(rand()%SIZE,rand()%SIZE);
    plates.push_back(plate);

  }

  //Generate Plate Centroids
  sample::disc(centroids, K, glm::vec2(0), glm::vec2(256));

  for(int i = 0; i < centroids.size(); i++){
    Litho newseg(0.5f, 0.0, &centroids[i]); //Properly Scaled Position
    segments.push_back(newseg);

    int nearest = rand()%nplates;
    float dist = SIZE*SIZE;
    for(int j = 0; j < nplates; j++){

      if(glm::length(plates[j].pos-centroids[i]) < dist){
        dist = glm::length(plates[j].pos-centroids[i]);
        nearest = j;
      }

    }

    //Add References
    plates[nearest].segments.push_back(i);

  }

  for(int j = 0; j < nplates; j++)
    plates[j].recenter(segments);

}

void World::cluster(Shader* voronoi, Instance* inst){

  clustering->target(glm::vec3(1));
  voronoi->use();
  voronoi->uniform("R", R);
  voronoi->uniform("depthmap", false);
  inst->render();

  depthmap->target(glm::vec3(1));
  voronoi->use();
  voronoi->uniform("R", R);
  voronoi->uniform("depthmap", true);
  inst->render();

}

void World::drift(Instance* inst){

  //This could be causing garbage collection issues, not sure
  int* c = clustering->sample<int>(glm::vec2(0), dim, GL_COLOR_ATTACHMENT0, GL_RGBA);

  for(auto& p: plates){
    p.convect(heatmap, segments);
    p.grow(heatmap, segments);
    p.collide(c, segments);
  }

  inst->updateBuffer(centroids, 0);

}

/*
    Construct a mesh from the world!
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

      if( aind < w->segments.size() )
        a = glm::vec3(i  , 10.0*w->segments[aind].height, j  );
      if( bind < w->segments.size() )
        b = glm::vec3(i  , 10.0*w->segments[bind].height, j+1);
      if( cind < w->segments.size() )
        c = glm::vec3(i+1, 10.0*w->segments[cind].height, j  );
      if( dind < w->segments.size() )
        d = glm::vec3(i+1, 10.0*w->segments[dind].height, j+1);

      glm::vec3 stonecolor = glm::vec3(0.8);
      glm::vec3 collidecolor = glm::vec3(0.0,1.0,0.0);
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
        if(w->segments[aind].colliding){
          m->colors.push_back(collidecolor.x);
          m->colors.push_back(collidecolor.y);
          m->colors.push_back(collidecolor.z);
        }
        else{
          m->colors.push_back(stonecolor.x);
          m->colors.push_back(stonecolor.y);
          m->colors.push_back(stonecolor.z);
        }
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      if(bind < w->segments.size()){
        if(w->segments[bind].colliding){
          m->colors.push_back(collidecolor.x);
          m->colors.push_back(collidecolor.y);
          m->colors.push_back(collidecolor.z);
        }
        else{
          m->colors.push_back(stonecolor.x);
          m->colors.push_back(stonecolor.y);
          m->colors.push_back(stonecolor.z);
        }        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      if(cind < w->segments.size()){
        if(w->segments[cind].colliding){
          m->colors.push_back(collidecolor.x);
          m->colors.push_back(collidecolor.y);
          m->colors.push_back(collidecolor.z);
        }
        else{
          m->colors.push_back(stonecolor.x);
          m->colors.push_back(stonecolor.y);
          m->colors.push_back(stonecolor.z);
        }        m->colors.push_back(1.0);
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
        if(w->segments[dind].colliding){
          m->colors.push_back(collidecolor.x);
          m->colors.push_back(collidecolor.y);
          m->colors.push_back(collidecolor.z);
        }
        else{
          m->colors.push_back(stonecolor.x);
          m->colors.push_back(stonecolor.y);
          m->colors.push_back(stonecolor.z);
        }
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      if(cind < w->segments.size()){
        if(w->segments[cind].colliding){
          m->colors.push_back(collidecolor.x);
          m->colors.push_back(collidecolor.y);
          m->colors.push_back(collidecolor.z);
        }
        else{
          m->colors.push_back(stonecolor.x);
          m->colors.push_back(stonecolor.y);
          m->colors.push_back(stonecolor.z);
        }
        m->colors.push_back(1.0);
      }
      else{
        m->colors.push_back(magmacolor.x);
        m->colors.push_back(magmacolor.y);
        m->colors.push_back(magmacolor.z);
        m->colors.push_back(1.0);
      }

      if(bind < w->segments.size()){
        if(w->segments[bind].colliding){
          m->colors.push_back(collidecolor.x);
          m->colors.push_back(collidecolor.y);
          m->colors.push_back(collidecolor.z);
        }
        else{
          m->colors.push_back(stonecolor.x);
          m->colors.push_back(stonecolor.y);
          m->colors.push_back(stonecolor.z);
        }
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
};

std::function<glm::vec4(double)> heatmap = [](double hm){
  return glm::mix(glm::vec4(0.0,0.0,0.0,1.0), glm::vec4(1.0), hm);
};
