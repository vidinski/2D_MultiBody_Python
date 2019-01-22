import numpy as np
from scipy.integrate import odeint
from body_coordinates import body_coordinates, BC_constraints, joints

g = 9.81
PI = np.pi

def solveSys(x,t,invM, K, bodies, tnb, joint_list):
    #_______________________________________________________________________
    x = np.matrix(x)
    #print(x[0,3*1:3*1+3])
    q = np.zeros([tnb,6]); 
    bodies[0].BC_trans(0.0, np.transpose([0.0, 0.0])) 

    for i in range (1,len(bodies)):
        j = i-1
        q[j,:] = np.concatenate((x[0,3*j:3*j+3],x[0,3*(tnb+j):3*(tnb+j)+3]), axis = 1)
        bodies[i].BC_trans(q[j,2], np.transpose([[q[j,0], q[j,1]]]))

    #Constraints Jacobian: 
    #_______________________________________________________________________
    [Jacij, gammaij] = BC_constraints(bodies[joint_list[0].body_i[0]],
		          	      joint_list[0].body_i[1],
		   		      bodies[joint_list[0].body_j[0]], 
		   		      joint_list[0].body_j[1],
		   		      joint_list[0].joint_type, 
				      np.matrix([0,0,0,0,0,0]),
                                      np.matrix(q[0,:]), 
    				      tnb)
    D = Jacij
    gamma = np.transpose(np.matrix(gammaij))

    if len(joint_list) > 1: 
        for i in range (1,len(joint_list)):  
	    #print(joint_list[i].body_i[1])
            [Jacij, gammaij] = BC_constraints(bodies[joint_list[i].body_i[0]],
		          	              joint_list[i].body_i[1],
		   		              bodies[joint_list[i].body_j[0]], 
		   		              joint_list[i].body_j[1],
		   		              joint_list[i].joint_type, 
		   		              np.matrix(q[joint_list[i].body_i[0]-1,:]), 
                                              np.matrix(q[joint_list[i].body_j[0]-1,:]),
                                              tnb)
	    D = np.concatenate((D,Jacij),axis = 0)
	    gammaij = np.transpose(np.matrix(gammaij))
	    gamma = np.concatenate((gamma,gammaij),axis = 0)
    F = forces(x,t,K,bodies)
    #F = np.transpose(F)
    DT = np.transpose(D)
    DMinv = np.matmul(D,invM)
    DMinvDT = np.matmul(DMinv,DT)
    DMinvF = np.matmul(DMinv, F)
    gamma_F = gamma-DMinvF
    lagr_lambda = np.linalg.solve(DMinvDT,gamma_F)
    qdoubledot = np.matmul(DT,lagr_lambda)
    qdoubledot = qdoubledot+F
    qdoubledot = np.matmul(invM, qdoubledot)   
    
    #extract velocities: 
    xyth_vel = x[0,tnb*3:tnb*6]
    qdoubledot = np.transpose(qdoubledot)
    xdot = np.concatenate((xyth_vel,qdoubledot),axis = 1 )
    xdot = np.array(xdot)
    return xdot[0]

def forces(x,t,K,bodies):
    g = -9.81*np.matrix([0.0,1.0,0.0])
    for i in range (1,len(bodies)):
        Fg = np.concatenate((g,g), axis = 1) 
    ### LQR CONTROL ###
    th = np.zeros(((len(bodies)-1)*2,1))
    th[0,0] =  (x[0,2]-PI*0.5)
    th[1,0] =  (x[0,5] - x[0,2])
    th[2,0] = x[0,8]    
    th[3,0] = x[0,11] - x[0,8]
    Fcontrol = [0.0, 0.0, np.matmul(K,th), 0.0, 0.0, -np.matmul(K,th)]
    #Step Force input
    if t < 0.1: 
       Fstep = [0.0, 0.0, 0.0, 0.30, 0.0, 0.0]
    #elif t >= 1.0 and t<1.01:
       #Fstep = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    else: 
       Fstep = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]    
    
    F = np.transpose(Fg+Fcontrol+Fstep)     
    #F = np.transpose(Fg+Fstep) 
    return F




