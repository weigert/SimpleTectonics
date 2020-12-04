#include "TinyEngine/TinyEngine.h"
#include "TinyEngine/include/helpers/color.h"
#include "TinyEngine/include/helpers/image.h"
#include "TinyEngine/include/helpers/timer.h"
#include <noise/noise.h>
#include <chrono>

#include "source/poisson.h"
#include "source/model.h"

#include "source/world.h"


int main( int argc, char* args[] ) {

	//Setup Window
	Tiny::view.vsync = false;
	Tiny::window("Plate Tectonics Simulation", WIDTH, HEIGHT);

	Tiny::event.handler  = eventHandler;
	Tiny::view.interface = [](){};

	//Generate Seeded World
	int SEED = time(NULL);
	if(argc == 2)
		SEED = std::stoi(args[1]);
	World world(SEED);

	//Setup Shaders
	Shader shader({"source/shader/default.vs", "source/shader/default.fs"}, {"in_Position", "in_Normal", "in_Color"});
	Shader depth({"source/shader/depth.vs", "source/shader/depth.fs"}, {"in_Position"});
	Shader billboardshader({"source/shader/flat.vs", "source/shader/flat.fs"}, {"in_Quad", "in_Tex"});

	Shader voronoi({"source/shader/voronoi.vs", "source/shader/voronoi.fs"}, {"in_Quad", "in_Tex", "in_Centroid"});
	Shader diffusion({"source/shader/flat.vs", "source/shader/diffusion.fs"}, {"in_Quad", "in_Tex"});
	Shader subduction({"source/shader/flat.vs", "source/shader/subduction.fs"}, {"in_Quad", "in_Tex"}, {"colliding"});
	Shader sedimentation({"source/shader/flat.vs", "source/shader/sedimentation.fs"}, {"in_Quad", "in_Tex"}, {"colliding"});

	//Utility Classes
	Square2D flat;
	Model model(tectonicmesh, &world);
	Billboard shadow(2000, 2000, true);			//Shadow Map

	//Prepare instance render of flat, per-centroid
	Instance instance(&flat);
	instance.addBuffer(world.centroids);

	int n = 0;
	float us = 0.0; //Rolling average execution time calculation in microseconds (us)

	world.cluster(&voronoi, &instance);
	model.construct(tectonicmesh, &world); //Reconstruct Updated Model

	Tiny::view.pipeline = [&](){

		//Render Shadowmap
		shadow.target();                  //Prepare Target
		depth.use();                      //Prepare Shader
		model.model = glm::translate(glm::mat4(1.0), -viewPos);
		depth.uniform("dmvp", depthProjection * depthCamera * model.model);
		model.render(GL_TRIANGLES);       //Render Model

		//Regular Image
		if(viewplates) Tiny::view.target(skyCol);
		else Tiny::view.target(skyBlue);

		shader.use();                   //Prepare Shader
		shader.texture("shadowMap", shadow.depth);
    shader.uniform("lightCol", lightCol);
    shader.uniform("lightPos", lightPos);
    shader.uniform("lookDir", lookPos-cameraPos);
    shader.uniform("lightStrength", lightStrength);
    shader.uniform("projectionCamera", projection * camera);
    shader.uniform("dbmvp", biasMatrix * depthProjection * depthCamera * glm::mat4(1.0f));
    shader.uniform("model", model.model);
    model.render(GL_TRIANGLES);    //Render Model

		if(viewmap){

			billboardshader.use();
			billboardshader.texture("imageTexture", world.clustering->texture);
			flat.move(glm::vec3(-1.0+0.25/WIDTH*HEIGHT,1.0-0.25,0.0), 0, glm::vec3(1.0f*0.25/WIDTH*HEIGHT,0.25,0.0));
			billboardshader.uniform("model", flat.model);
			flat.render();

			billboardshader.use();
			billboardshader.texture("imageTexture", world.heatA->texture);
			flat.move(glm::vec3(-1.0+0.75/WIDTH*HEIGHT,1.0-0.25,0.0), 0, glm::vec3(1.0f*0.25/WIDTH*HEIGHT,0.25,0.0));
			billboardshader.uniform("model", flat.model);
			flat.render();

			billboardshader.use();
			billboardshader.texture("imageTexture", world.heightA->texture);
			flat.move(glm::vec3(-1.0+1.25/WIDTH*HEIGHT,1.0-0.25,0.0), 0, glm::vec3(1.0f*0.25/WIDTH*HEIGHT,0.25,0.0));
			billboardshader.uniform("model", flat.model);
			flat.render();

		}

	};

	int m = 0;

	Tiny::loop([&](){ //Execute every frame

		if(animate){

			world.drift();
			world.diffuse(&diffusion, &subduction, &flat);
			world.addRock(&diffusion, &sedimentation, &flat);
			world.update(&instance);
			world.cluster(&voronoi, &instance);
			model.construct(tectonicmesh, &world); //Reconstruct Updated Model

		}

	});

	Tiny::quit();

	return 0;
}
