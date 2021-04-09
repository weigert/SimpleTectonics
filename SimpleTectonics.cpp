#include <TinyEngine/TinyEngine>
#include <TinyEngine/color>
#include <TinyEngine/image>
#include <TinyEngine/timer>
#include <noise/noise.h>
#include <chrono>

//#define K 16384
#define K 8192
//#define K 4096
//#define K 2048

#define SIZE 256
#define nplates 12
#define DT 0.025f
#define WIDTH 1000
#define HEIGHT 1000
const float R = 2.0f*sqrt(4.0f/3.14159265f/K);
const float CR = R/14.0f;	//Collision Radius

#include "source/poisson.h"
#include "source/cluster.h"
#include "source/tectonics.h"
#include "source/scene.h"

int main( int argc, char* args[] ) {

	Tiny::view.vsync = false;
	Tiny::window("Plate Tectonics Simulation", WIDTH, HEIGHT);

	Tiny::event.handler  = eventHandler;
	Tiny::view.interface = interfaceFunc;

	int SEED = time(NULL);
	if(argc == 2)
		SEED = std::stoi(args[1]);

	//Setup Shaders
	Shader shader({"source/shader/default.vs", "source/shader/default.fs"}, {"in_Position", "in_Normal", "in_Color"});
	Shader depth({"source/shader/depth.vs", "source/shader/depth.fs"}, {"in_Position"});
	Shader billboardshader({"source/shader/flat.vs", "source/shader/flat.fs"}, {"in_Quad", "in_Tex"});

	Shader diffusion({"source/shader/flat.vs",  "source/shader/diffusion.fs"},  {"in_Quad", "in_Tex"});
	Shader cascading({"source/shader/flat.vs",  "source/shader/cascading.fs"},  {"in_Quad", "in_Tex"}, {"height"});
	Shader subduction({"source/shader/flat.vs", "source/shader/subduction.fs"}, {"in_Quad", "in_Tex"}, {"colliding"});

	World world(SEED);
	world.heightB->target(vec3(0.3));	//Sea-Level
	world.heightA->target(vec3(0.3));

	//Utility Classes
	Square2D flat;
	Model platemodel(tectonicmesh, &world);
	Model earthmodel(tectonicmesh, &world);
	Billboard shadow(2000, 2000, true);			//Shadow Map

/*
	int n = 500;
	while(n-- > 0){

		//Move the Plates
		for(auto& p: world.plates)
			p.update(world.cluster, world.heatmap);

		//Compute Effect on Heat and Mass
		world.subduct(&diffusion, &subduction, &flat, 25);

	//	if(viewsurface)

		//Update the World
		world.update();
		std::cout<<"CYCLE "<<n<<std::endl;

	}

	viewsurface = true;

world.sediment(&cascading, &flat, 25);
viewplates = false;
earthmodel.construct(tectonicmesh, &world); //Reconstruct Updated Model

	animate = false;

*/


/*

*/

	//Construct the Models
	viewplates = true;
	platemodel.construct(tectonicmesh, &world);

	if(viewsurface){
		viewplates = false;
		earthmodel.construct(tectonicmesh, &world);
	}

	Tiny::view.pipeline = [&](){

		if(viewsurface){

			shadow.target(true);                  																//Prepare Target
			depth.use();                      																		//Prepare Shader
			earthmodel.model = glm::translate(glm::mat4(1.0), -viewPos + vec3(0,0,0));
			depth.uniform("dmvp", depthProjection*depthCamera*earthmodel.model);
			earthmodel.render(GL_TRIANGLES);       																		//Render Model

			Tiny::view.target(skyBlue);

			shader.use();                  																				//Prepare Shader
			shader.texture("shadowMap", shadow.depth);
			shader.uniform("lightCol", lightCol);
			shader.uniform("lightPos", lightPos);
			shader.uniform("lookDir", lookPos-cameraPos);
			shader.uniform("lightStrength", lightStrength);
			shader.uniform("projectionCamera", projection * camera);
			shader.uniform("dbmvp", biasMatrix*depthProjection*depthCamera);
			shader.uniform("model", earthmodel.model);
			earthmodel.render(GL_TRIANGLES);

		}

		if(viewplates){
			shadow.target(true);                  																		//Prepare Target
			depth.use();                      																		//Prepare Shader
			platemodel.model = glm::translate(glm::mat4(1.0), -viewPos + vec3(0,-75,0));
			depth.uniform("dmvp", depthProjection*depthCamera*platemodel.model);
			platemodel.render(GL_TRIANGLES);       																		//Render Model

			Tiny::view.target(skyBlue);

			shader.use();                  																				//Prepare Shader
			shader.texture("shadowMap", shadow.depth);
			shader.uniform("lightCol", lightCol);
			shader.uniform("lightPos", lightPos);
			shader.uniform("lookDir", lookPos-cameraPos);
			shader.uniform("lightStrength", lightStrength);
			shader.uniform("projectionCamera", projection * camera);
			shader.uniform("dbmvp", biasMatrix*depthProjection*depthCamera);
			shader.uniform("model", platemodel.model);
			platemodel.render(GL_TRIANGLES);    																				//Render Model
		}

		if(viewmap){

			billboardshader.use();
			billboardshader.texture("imageTexture", world.cluster.target->texture);
			flat.move(glm::vec3(-1.0+0.25/WIDTH*HEIGHT,1.0-0.25,0.0), 0, glm::vec3(1.0f*0.25/WIDTH*HEIGHT,0.25,0.0));
			billboardshader.uniform("model", flat.model);
			flat.render();

			flat.model = glm::mat4(1);


			billboardshader.use();
			billboardshader.texture("imageTexture", world.heatA->texture);
			flat.move(glm::vec3(-1.0+0.75/WIDTH*HEIGHT,1.0-0.25,0.0), 0, glm::vec3(1.0f*0.25/WIDTH*HEIGHT,0.25,0.0));
			billboardshader.uniform("model", flat.model);
			flat.render();

			flat.model = glm::mat4(1);

			billboardshader.use();
			billboardshader.texture("imageTexture", world.heightA->texture);
			flat.move(glm::vec3(-1.0+1.25/WIDTH*HEIGHT,1.0-0.25,0.0), 0, glm::vec3(1.0f*0.25/WIDTH*HEIGHT,0.25,0.0));
			billboardshader.uniform("model", flat.model);
			flat.render();

		}

	};

	Tiny::loop([&](){ //Execute every frame

		if(animate){

			//Move the Plates
		  for(auto& p: world.plates)
				p.update(world.cluster, world.heatmap);

			//Compute Effect on Heat and Mass
			world.subduct(&diffusion, &subduction, &flat, 25);

			if(viewsurface) world.sediment(&cascading, &flat, 25);

			//Update the World
			world.update();

			//Construct the Models
			if(!viewsurface){
				viewplates = true;
				platemodel.construct(tectonicmesh, &world);
			}

			if(viewsurface){
				viewplates = false;
				earthmodel.construct(tectonicmesh, &world);
			}

		}

	});

	Tiny::quit();

	return 0;
}
