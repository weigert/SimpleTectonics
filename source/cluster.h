/*
================================================================================
                            Clustered Convection
================================================================================
*/

#include <vector>
using namespace glm;
using namespace std;

struct Segment {
  Segment(vec2* p):pos{p}{}
  vec2* pos = NULL;           //Segment Position
  int area = 1.0;             //Area of Segment (Pixels)
};

template<typename T>
class Cluster {
public:

Cluster(){
  init();
}

vector<vec2> points;     //Raw Centroid Data
vector<T*> segs;

int* indexmap;

Square2D* flat;
Instance* instance;

Billboard* target;
Shader* voronoi;

void update(){

  //SSBO for Area Accumulation
  vector<int> area;
  for(size_t i = 0; i < segs.size(); i++)
    area.push_back(0);
  voronoi->buffer("area", area);

  //Update Instance with Point-Set
  instance->updateBuffer(points, 0);

  //Render Voronoi Texture
  target->target(vec3(1));
  voronoi->use();
  voronoi->uniform("R", R);
  voronoi->uniform("depthmap", false);
  instance->render();

  //Extract the Index Map
  target->sample<int>(indexmap, vec2(0), vec2(SIZE, SIZE),
                      GL_COLOR_ATTACHMENT0, GL_RGBA);

  //Extract per-segment Area
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, voronoi->ssbo["area"]);
  int* areapointer = (int*)glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY);
  for(size_t i = 0; i < segs.size(); i++)
    segs[i]->area = segs[i]->area*0.99 + 0.01*areapointer[i];
  glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

}

T* add(vec2& p){

  T* newseg = new T(&p); //Properly Scaled Position
  segs.push_back(newseg);
  return newseg;

}

void reassign(){

  for(int i = 0; i < segs.size(); i++)
    segs[i]->pos = &points[i];

}

//Remove Matching Segments...
void remove(std::function<bool(T*)> f){

  for(int i = 0; i < segs.size(); i++){
    if(!f(segs[i])) continue;
    points.erase(points.begin()+i);
    delete segs[i];
    segs.erase(segs.begin()+i--);
  }
  reassign();

}

void init(){

  flat = new Square2D();                  //Clustering Render Model
  instance = new Instance(&*flat);        //Instance of Render Model
  target = new Billboard(SIZE, SIZE);     //Clustering Render Target
  indexmap = new int[SIZE*SIZE];          //Final Extracted Indexmap
  voronoi = new Shader( {"source/shader/voronoi.vs", "source/shader/voronoi.fs"},
                        {"in_Quad", "in_Tex", "in_Centroid"}, {"area"});

  sample::disc(points, K, vec2(0), vec2(256));
  for(auto&c : points) add(c);
  instance->addBuffer(points);

  update();

}

void quit(){

  delete target;
  delete[] indexmap;
  delete voronoi;

  for(int i = 0; i < segs.size(); i++)
    delete segs[i];

}

//Get the Segment Index at position p
int sample(vec2 p){

  int cmind = (int)p.y*SIZE+(int)p.x;
  vec4 col = color::i2rgba(indexmap[cmind]);  //Extract the Color
  if(col.z == 255) return -1;
  return (int)col.x + (int)col.y*256 + (int)col.z*256*256;

}

};
