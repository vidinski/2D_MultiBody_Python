#Phillip Vidinski 
#January 22nd 2019
# 2D_Multibody_Python

This version only supports revolute joints. 
Handling of forces can be ameliorated 

Script uses the body coordinate method and was written using the contents from: 

"Planar Multibody Dynamics: Forumlation, Programming, and Application" 
by Parviz E. Nikravesh


******************************************************************************************
The current code presents the acrobot problem 
(double pendulum actuated at the "elbow" joint). 
******************************************************************************************
Edit Main.py for setting up bodies and joints. Edit the forces in the "forces" function of solver.py
The animation can then be ran. 

*******************************************************************************************
Edit the links' center of mass, initial conditions, etc. example: 
*******************************************************************************************

*******************************************************************************************
#SETUP BODIES
*******************************************************************************************

####GROUND#####
body0 = body_coordinates(0, #index
			 np.transpose([[0.0, 0.0]]), #center
			 np.transpose([[0.0, 0.0]]), #shape
			 np.transpose([[0.0, 0.0],[0.0,0.0]]), #joints
			 np.transpose([[0.0, 1.0], [1.0,0.0]]), #unit vectors
			 0.0, #angle
			 100000.0, #mass
			 100000.0) #inertia 
bodies.append(body0)
		
##_______________________________________________________________________________________##
body1 = body_coordinates(1, #index
			 np.transpose([[0.0,0.5]]), #center
			 np.transpose([[0.5,-0.1],[-0.5,-0.1],[-0.5, 0.1],[0.5,0.1]]), #shape
			 np.transpose([[-0.5,0.0],[0.5,0.0]]), #joints
			 np.transpose([[0.0,1.0],[1.0, 0.0]]), #unit vectors
			 0.5*PI, #angle
			 1.0, #mass
			 1.0/12.0) #inertia
bodies.append(body1)
go.append(patches.Polygon(np.array(np.transpose(bodies[1].xy_local_shape))))


*******************************************************************************************
Edit the joits example: 
*******************************************************************************************

*******************************************************************************************
#SETUP JOINTS
*******************************************************************************************


joint_list = []

joint1 = joints("revolute", 1, [body0.index, 0], [body1.index,0])

joint_list.append(joint1)





