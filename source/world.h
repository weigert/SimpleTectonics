
//Plate Centroid
struct Litho{
public:

  Litho(glm::vec2 p){
    pos = p;
  }

  Litho(float d, float t, glm::vec2 p){
    density = d; thickness = t; pos = p;
  }

  float density = 1.0f;
  float thickness = 0.0f;
  float height = 0.0f;
  glm::vec2 pos;

};

//Prepare Noise for jiggling the centroids

class World {
public:

  World(){

    SEED = time(NULL);
    srand(SEED);

    perlin.SetOctaveCount(4);
    perlin.SetFrequency(1.0);
    perlin.SetPersistence(0.5);

    initialize();

    clustering = new Billboard(SIZE, SIZE);

  }

  ~World(){
    delete clustering;
  }

  //General Information
  int SEED = 0;
  float t = 0.0;
  const glm::vec2 dim = glm::vec2(SIZE, SIZE);
  noise::module::Perlin perlin;

  //Plate Centroids
  std::vector<glm::vec2> centroids; //Raw Position Buffer
  std::vector<glm::vec2> offset; //Raw Position Buffer
  std::vector<Litho> plates;        //Additional Data
  Billboard* clustering;
  void initialize();
  void drift(Instance* inst);
  void cluster(Shader* voronoi, Instance* inst);

};


void World::initialize(){

  //Generate Plate Centroids
  sample::disc(centroids, K, glm::vec2(-1), glm::vec2(1));
  offset = centroids;

  for(auto&c: centroids){
    float t = (float)(rand()%100)/100.0f;
    Litho newplate(1.0f, t, c);
    plates.push_back(newplate);
  }

  //Basically perform the clustering in the texture, then use that to construct the world.
  //Then choose the height as the centroid index
}

void World::cluster(Shader* voronoi, Instance* inst){

  clustering->target(color::black);
  voronoi->use();
  voronoi->uniform("R", R);
  inst->render();

}

void World::drift(Instance* inst){

  t += 0.005;

	for(unsigned int i = 0; i < centroids.size(); i++){
		offset[i].x = centroids[i].x + 0.5f*R*perlin.GetValue(centroids[i].x, centroids[i].y, t);
		offset[i].y = centroids[i].y + 0.5f*R*perlin.GetValue(centroids[i].x, centroids[i].y, -t);
	}

  inst->updateBuffer(offset, 0);

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
      glm::vec3 a = glm::vec3(i  , 5.0*w->plates[aind].thickness, j  );
      glm::vec3 b = glm::vec3(i  , 5.0*w->plates[bind].thickness, j+1);
      glm::vec3 c = glm::vec3(i+1, 5.0*w->plates[cind].thickness, j  );
      glm::vec3 d = glm::vec3(i+1, 5.0*w->plates[dind].thickness, j+1);

      //UPPER TRIANGLE

      //Get the Color of the Ground (Water vs. Flat)
      glm::vec3 color = glm::vec3(0.8);

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

      glm::vec3 n1 = glm::normalize(glm::cross(a-b, c-b));

      for(int i = 0; i < 3; i++){
        m->normals.push_back(n1.x);
        m->normals.push_back(n1.y);
        m->normals.push_back(n1.z);
        m->colors.push_back(color.x);
        m->colors.push_back(color.y);
        m->colors.push_back(color.z);
        m->colors.push_back(1.0);
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

      glm::vec3 n2 = glm::normalize(glm::cross(d-c, b-c));

      for(int i = 0; i < 3; i++){
        m->normals.push_back(n2.x);
        m->normals.push_back(n2.y);
        m->normals.push_back(n2.z);
        m->colors.push_back(color.x);
        m->colors.push_back(color.y);
        m->colors.push_back(color.z);
        m->colors.push_back(1.0);
      }

    }
  }
};
