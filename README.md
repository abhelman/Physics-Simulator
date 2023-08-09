# Physics-Simulator

Physics Simulator built using openGL.

# Walking Robot
A 4-legged robot built out of springs and masses simulates walking. 

The colors of the robot represent different materials. Dark blue is bone, light blue is tissue, and pink/green are contracting/expanding muscle. For each material, a different spring constant is used, with bone being very rigid, muscles medium, and tissue stretchy. The "breathing" of the springs is given by the equation: 
l = l_init*(1.0 + b*sin (w*t + c)) where w = 6.28, and (b, c, k) are as follows: 1. Tissue: (0, 0, 1000) 2. Muscle1: (.25, 0, 5000) 3. Muscle2: (.25, 3.14, 5000) 4. Bone: (0, 0, 20000). 

An evolutionary algorithm was used to determine the material for each spring.

https://github.com/abhelman/Physics-Simulator/assets/113809366/1fbbf04c-8c5d-49da-8c77-b8dee65bf501

# Bouncing Cube
A rigid cube given a small spin bounces on a hard surface.

https://github.com/abhelman/Physics-Simulator/assets/113809366/1e51e3b1-9eba-4f57-81c4-697abc56857f

# Sources

All graphics were created using the OpenGl graphics library and the learnopengl.com tutorial. The following API were used to implement OpengGL:
	GLFW, v.3.3.5, 28/10/21; https://www.glfw.org/download.html
 	API, Glad API, and OpenGL Mathematics (GML); https://github.com/g-truc/glm
	Glad, v.3.3; https://glad.dav1d.de/
 
To create the spheres and rods, the below source was used as a reference:
		Song Ho Ah, 01/11/17; http://www.songho.ca/opengl/gl_sphere.html
  
The source code for the Titan library was consulted to determine correct use of friction:
	J. Austin, R. Corrales-Fatou, S. Wyetzner, and H. Lipson, “Titan: A Parallel Asynchronous Library for Multi-Agent and Soft-Body Robotics using NVIDIA CUDA,” ICRA 2020, May 2020, https://github.com/jacobaustin123/Titan/blob/master/src/object.cu
