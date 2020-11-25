class Target {
public:
  Target(){ glGenFramebuffers(1, &fbo); };              //Default constructor

  Target(int W, int H, bool c = false, bool d = true):  //Construct with a size
    Target(){ WIDTH = W; HEIGHT = H; dAttach = d; cAttach = c; }

  ~Target(){ glDeleteFramebuffers(1, &fbo); }           //Default destructor

  unsigned int WIDTH, HEIGHT;
  bool dAttach = false, cAttach = true;                 //Whether we have a depth / color attachment
  GLuint fbo;                                           //FBO (OpenGL: everything is an int)

  template<typename T>
  void bind(T& t, bool d){                              //Bind a texture to the FBO
    if(d) glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, t.texture, 0);
    else  glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, t.texture, 0);
  }

  template<typename T>
  void setup(T& t, T& d){                               //Add color and depth textures!
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    if(dAttach) bind(d, true);
    if(cAttach) bind(t, false);
    else{ //If we don't have a color attachment, don't draw to it
      glDrawBuffer(GL_NONE);
      glReadBuffer(GL_NONE);
    }

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      std::cout<<"Framebuffer Incomplete"<<std::endl;

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
  }

  void target(bool t = true){          //This FBO as a render target (no clearing)
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0, 0, WIDTH, HEIGHT);
    if(t && dAttach) glClear( GL_DEPTH_BUFFER_BIT );
    if(t && cAttach) glClear( GL_COLOR_BUFFER_BIT );
  }
  template<typename T> void target(T a){
    glClearColor(a[0], a[1], a[2], 1.0f);
    target();
  }

  template<typename T> //Raw Buffer Sampling
  T* raw(glm::vec2 p, glm::vec2 d, GLenum ATTACH, GLenum FORMAT){
    if(d.x <= 0 || d.y <= 0 || p.x+d.x > WIDTH || p.y+d.y > HEIGHT) return NULL;
    T* h = new T[(int)(d.x*d.y)];
    glBindFramebuffer(GL_FRAMEBUFFER, fbo); glReadBuffer(ATTACH);
    glReadPixels(p.x, p.y, d.x, d.y, FORMAT, GL_UNSIGNED_BYTE, h);
    return h;
  }
  char* sampledepth(glm::vec2 p, glm::vec2 d){
    return raw<char>(p, d, GL_DEPTH_ATTACHMENT, GL_DEPTH_COMPONENT);
  }
  int* sample(glm::vec2 p, glm::vec2 d){
    return raw<int>(p, d, GL_COLOR_ATTACHMENT0, GL_RGBA);
  }
  glm::vec4 sample(glm::vec2 p){
    int* val = sample(p, glm::vec2(1));
    float r = ((*val >>  0) & 0xff);
    float g = ((*val >>  8) & 0xff);
    float b = ((*val >> 16) & 0xff);
    float a = ((*val >> 24) & 0xff);
    return glm::vec4(r,g,b,a);
  }

};

class Billboard: public Target{   //Billboard specialization
public:
  Texture texture, depth;         //Two normal textures

  Billboard(int W, int H, bool c = true, bool d = true):
  Target(W, H, c, d){
    if(dAttach) depth.depth(WIDTH, HEIGHT);
    if(cAttach) texture.empty(WIDTH, HEIGHT);
    setup(texture, depth);        //Bind the two normal textures to the billboard
  }

  Billboard(SDL_Surface* s):      //Render target preloaded with an image
  Billboard(s->w, s->h, true, false){
    texture.raw(s);               //Construct the texture from raw surface data
    bind(texture, false);         //Bind it to the billboard as color texture
  }
};

class Cubemap: public Target{     //Cubemap specialization
public:
  Cubetexture texture, depth;     //Two cubemap textures

  Cubemap(int W, int H, bool c = true, bool d = true):Target(W, H, c, d){
    if(dAttach) depth.depth(WIDTH, HEIGHT);
    if(cAttach) texture.empty(WIDTH, HEIGHT);
    setup(texture, depth);
  }
};
