#include "TinyEngine/TinyEngine.h"
#include "TinyEngine/include/helpers/color.h"
#include "TinyEngine/include/helpers/image.h"
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
	Shader voronoi({"source/shader/voronoi.vs", "source/shader/voronoi.fs"}, {"in_Quad", "in_Tex", "in_Centroid"});
	Shader billboardshader({"source/shader/billboard.vs", "source/shader/billboard.fs"}, {"in_Quad", "in_Tex"});
	Shader shader({"source/shader/default.vs", "source/shader/default.fs"}, {"in_Position", "in_Normal", "in_Color"});
	Shader depth({"source/shader/depth.vs", "source/shader/depth.fs"}, {"in_Position"});

	//Utility Classes
	Square2D flat;
	Model model(tectonicmesh, &world);

	Billboard shadow(2000, 2000, true);
	Billboard image(WIDTH, HEIGHT, false); //1200x800, depth only

	//Prepare instance render of flat, per-centroid
	Instance instance(&flat);
	instance.addBuffer(world.centroids);

	int n = 0;
	float us = 0.0; //Rolling average execution time calculation in microseconds (us)

	world.cluster(&voronoi, &instance);
	model.construct(tectonicmesh, &world); //Reconstruct Updated Model

	//Setup 2D Images
  Billboard map(image::make<double>(world.dim, world.heatmap, heatmap));

	Tiny::view.pipeline = [&](){

		/*
				Convert this into a heightmap!
		*/

		//Render Shadowmap
		shadow.target();                  //Prepare Target
		depth.use();                      //Prepare Shader
		model.model = glm::translate(glm::mat4(1.0), -viewPos);
		depth.uniform("dmvp", depthProjection * depthCamera * model.model);
		model.render(GL_TRIANGLES);       //Render Model

		//Regular Image
    //image.target(skyCol);           //Prepare Target
		Tiny::view.target(skyCol);
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

/*
    //Render to Screen
    Tiny::view.target(color::black);    //Prepare Target
    effect.use();                //Prepare Shader
		effect.texture("imageTexture", image.texture);
		effect.texture("depthTexture", image.depth);
		effect.uniform("model", flat.model);
		flat.render();
    //image.render();                     //Render Image
*/

		if(viewmap){

			billboardshader.use();

			billboardshader.texture("imageTexture", world.clustering->texture);
			flat.move(glm::vec3(-1.0+0.25/WIDTH*HEIGHT,1.0-0.25,0.0), 0, glm::vec3(1.0f*0.25/WIDTH*HEIGHT,0.25,0.0));
			billboardshader.uniform("model", flat.model);
			flat.render();

			billboardshader.texture("imageTexture", map.texture);
			flat.move(glm::vec3(-1.0+0.75/WIDTH*HEIGHT,1.0-0.25,0.0), 0, glm::vec3(1.0f*0.25/WIDTH*HEIGHT,0.25,0.0));
			billboardshader.uniform("model", flat.model);
			flat.render();

		}


	};

	int m = 0;

	Tiny::loop([&](){ //Execute every frame


		if(animate){

			world.drift(&instance);
			world.cluster(&voronoi, &instance);
			model.construct(tectonicmesh, &world); //Reconstruct Updated Model

		}


	});

	Tiny::quit();

	return 0;
}
