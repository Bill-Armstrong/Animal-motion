The idea of this project is to challenge the machine-learning community to do as well as a new-born animal. It can quickly learn to control its limbs, get up off the ground, move around and then run like the wind. 

This **C++** program was developed using **QtCreator** and **Qt3D** and also employs Eigen matrices and quaternions for the dynamics. It is based upon a recursive method for solving the equations of motion of the robot arm on the space shuttle, later generalized to tree-linkages and applied to computer animation as described in the paper by William W. Armstrong and Mark W. Green, The dynamics of articulated rigid bodies for purposes of animation, Visual Computer (1985) 1(231:240) Springer-Verlag.

The current status is that the simulation works and shows five cylinders forming a tree-linkage supposed to look like a four-footed animal.  The later addition of links for head, neck and other extremities requires very few changes to the program, but the control will be a lot harder, so it is not done here.

The computer graphics need some work. I tried to make a movie, but all I have been able to do is get one chosen frame to be shown.  The movie is really necessary to be able to judge whether the motion is realistic.

When the above is working, I intend to add a machine learning component based on the ALN-machine-learning repository on Git Hub.  Control is currently done by imposing desired orientations of the four links representing legs, with a restorative torque. There is no provision for active control, but the "animal" can fall on the ground and bounce around a bit.

The ultimate goal is to support research into android robot motion or prostheses that help persons with spinal cord injuries or other impediments to motion. Collaborators are welcome.

Bill Armstrong
