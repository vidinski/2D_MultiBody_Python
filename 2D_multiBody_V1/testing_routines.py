"""
############################################################################################
#TEST SOLVER ROUTINE
q = np.matrix([1.0,2.0,3.0,4.0,5.0,6.0]) 
print(q[0,:][0,3:6])
#q = np.zeros([tnb,6]); 
#print(q[0,:][0,5])
#print(q)

#q = np.transpose(q)
bodies[0].BC_trans(0.0, np.transpose([0.0, 0.0])) 

for i in range (1,len(bodies)):
    j = i-1
    q[j,:] = np.concatenate((x0[3*j:3*j+3],x0[3*(tnb+j):3*(tnb+j)+3]), axis = 0)
    bodies[i].BC_trans(q[j,2], np.transpose([[q[j,0], q[j,1]]]))

#Constraints Jacobian: 
#_______________________________________________________________________
[Jacij, gammaij] = BC_constraints(bodies[joint_list[0].body_i[0]],
		          	  joint_list[0].body_i[1],
		   		  bodies[joint_list[0].body_j[0]], 
		   		  joint_list[0].body_j[1],
		   		  joint_list[0].joint_type, 
				  q[0,:],
				  np.matrix(np.zeros([1,6])),
				  tnb)
D = Jacij
gamma = np.transpose(np.matrix(gammaij))

if len(joint_list) > 1: 
    for i in range (1,len(joint_list)):    
        [Jacij, gammaij] = BC_constraints(bodies[joint_list[i].body_i[0]],
		          	          joint_list[i].body_i[1],
		   		          bodies[joint_list[i].body_j[0]], 
		   		          joint_list[i].body_j[1],
		   		          joint_list[i].joint_type, 
		   		          q[joint_list[i].body_i[0]], 
                                          q[joint_list[i].body_j[0],:], 
                                          tnb)
	D = np.concatenate((D,Jacij),axis = 0)
	gammaij = np.transpose(np.matrix(gammaij))
	gamma = np.concatenate((gamma,gammaij),axis = 0)
DT = np.transpose(D)
DMinvDT = np.matmul(D,invM)
DMinvDT = np.matmul(DMinvDT,np.transpose(D))
lagr_lambda = np.linalg.solve(DMinvDT,gamma)
qdoubledot = np.matmul(DT,lagr_lambda)
qdoubledot = np.matmul(invM,qdoubledot)




print(lagr_lambda)
print(gamma)
print(qdoubledot)
#Constraints Jacobian: 
#_______________________________________________________________________

"""






#"Test constraints" 
#________________________________________________________________________________
#q0 = [0.0, 0.0, 0.0,0.0, 0.0, 0.0]
#q1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#x1 = np.ones([20,20])
#q = np.zeros([2,5])
#q[0,:] = np.concatenate((x1[0,1:3],x1[1,1:4]), axis = 0)
#for i in range (0,(len(bodies)-1)):      
#    q[i,:] = np.concatenate((x1[i*3+2,1:3],x1[1,1:4]), axis = 0)

#print(q); 

#body1.BC_trans(qi[2], np.transpose([[qi[0], qi[1]]]))
#body2.BC_trans(qj[2], np.transpose([[qj[0], qj[1]]]))
#print(i)
#print(Jac)
#print(bodies[1].xy_global_center)
#[Jac1, gamma01] = BC_constraints(bodies[0], 0, bodies[1], 0, "revolute", q0, q1, tnb)
#print(body1.xy_global_center)
#print(body1.xy_jointrot)
#print(for_loop_bod.xy_global_center)
#print(Jac1)
#print(gamma01)
#print(joint1.joint_type)
