import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import matplotlib.patches as patches
import scipy.integrate as integrate
from solver import solveSys
from body_coordinates import body_coordinates, BC_constraints, joints
PI = np.pi








#go2 = patches.Polygon(np.transpose(body1.xy_global_shape),closed=True, fc='r', ec='r')
#patch3 =ax.add_patch(go2)

def init_anim(bodies):
    bodies[1].BC_trans(bodies[1].angle, bodies[1].xy_global_center)
    print(bodies[1].xy_global_shape+bodies[1].xy_global_center)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-3,13)
    ax.set_xlim(-10,10)

    background1 = patches.Rectangle([-50, 0], 100.0, 50.0, 0.0, color = 'blue', alpha = 0.35)
    background2 = patches.Rectangle([-50,-50], 100.0, 50.0, 0.0, color = 'green', alpha = 0.5)
    background1 = ax.add_patch(background1)
    background2 = ax.add_patch(background2)
    
    #for j in range (1,len(bodies)):
    
    return background1, background2, fig









