import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import scipy.integrate as integrate
from solver import solveSys
from body_coordinates import body_coordinates, BC_constraints, joints
PI = np.pi
bodies = []
go = []
tnb = 2
dt = 0.05 
t = np.arange(0.0, 10.0, dt)
#K = [[-6872.45297817, -2402.45660722, -2041.72935879,  -583.66541762]]

K = [[-6872.45297817*0.0, -2402.45660722*0.0, -2041.72935879,  -583.66541762]]
###########################################################################################
					#SETUP BODIES
###########################################################################################

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

##_______________________________________________________________________________________##
body2 = body_coordinates(2, #index
			 np.transpose([[0.0, 0.5*3.0]]), #center
			 np.transpose([[0.5,-0.1],[-0.5,-0.1],[-0.5, 0.1],[0.5,0.1]]), #shape
			 np.transpose([[-0.5,0.0],[0.5,0.0]]), #joints
			 np.transpose([[0.0,1.0],[1.0, 0.0]]), #unit vectors
			 0.5*PI,   #angle
			 1.0,      #mass
			 1.0/12.0) #inertia
bodies.append(body2)
go.append(patches.Polygon(np.array(np.transpose(bodies[2].xy_local_shape))))


###########################################################################################
					#SETUP ANIMATION WINDOW
###########################################################################################

fig = plt.figure()
#ax1 = fig.add_subplot(121)
ax1 = fig.add_subplot(111)
ax1.set_ylim(-3,4)
ax1.set_xlim(-3.5,3.5)
background1 = patches.Rectangle([-50, 0], 100.0, 50.0, 0.0, color = 'blue', alpha = 0.35)
background2 = patches.Rectangle([-50,-50], 100.0, 50.0, 0.0, color = 'green', alpha = 0.5)
background1 = ax1.add_patch(background1)
background2 = ax1.add_patch(background2)

for j in range (1,len(bodies)):
    jj = j-1
    bodies[j].BC_trans(bodies[j].angle,bodies[j].xy_global_center)
    go[jj].set_xy(np.transpose(bodies[j].xy_global_shape))
    patch = ax1.add_patch(go[jj])

###########################################################################################
					#SETUP JOINTS
###########################################################################################
joint_list = []

joint1 = joints("revolute", 1, [body0.index, 0], [body1.index,0])

joint_list.append(joint1)

##_______________________________________________________________________________________##
joint2 = joints("revolute", 2, [body1.index, 1], [body2.index,0])

joint_list.append(joint2)

###########################################################################################
			#SETUP EQUATIONS OF MOTION AND SOLVE 
###########################################################################################
M = bodies[1].mass
x0 = np.concatenate((bodies[1].xy_global_center, np.matrix(bodies[1].angle)), axis = 0)
if len(bodies) > 2:
    for i in range (2,len(bodies)):
        M = np.concatenate((M, bodies[i].mass), axis = 0)
        x0 = np.concatenate((x0,np.concatenate((bodies[i].xy_global_center, np.matrix(bodies[i].angle)), axis = 0) ), axis = 0)
M = np.diag(M)
invM = np.linalg.inv(M)
x0 = np.concatenate((np.transpose(x0),np.zeros([1,3*tnb])), axis = 1)
x0 = np.array(x0) 
#print(joint_list[1].body_j[0])
x = integrate.odeint(solveSys, x0[0],t, args = (invM,K,bodies,tnb,joint_list))

###########################################################################################
					#ANIMATION
###########################################################################################

def init():
    return [] #patch,
np.matrix(x)
def animate(i):
    q = np.zeros([tnb,6])   
    x_anim = x[i,:] 
    x_anim = np.matrix(x_anim)  
    for j in range (1,len(bodies)):
        jj = j-1
        q[jj,:] = np.concatenate((x_anim[0,3*jj:3*jj+3],x_anim[0,3*(tnb+jj):3*(tnb+jj)+3]),axis=1)
        bodies[j].BC_trans(q[jj,2], np.transpose([[q[jj,0], q[jj,1]]]))
        go[jj].set_xy(np.transpose(bodies[j].xy_global_shape))
        patch = ax1.add_patch(go[jj]) 
    return [] #patch,


ani = animation.FuncAnimation(fig, animate, np.arange(1, len(x)),
                              interval=1, blit=True, init_func=init)


#ax2 = fig.add_subplot(122)
#plt.plot(t,x[:,0:6])
plt.grid(True)
#plt.legend(('x','y','th','xd','yd','thd'))
plt.show()











