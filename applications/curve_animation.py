import sys
sys.path.append('./..')
from interfaces import *
from utilities import *

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from lib.progress_bar import *

############
# SETTINGS #
############

# for mac users it is required:
#   sudo port install ffmpeg

# this is the number of frames in the video
nframes = 200

################
# VECTOR SPACE #
################

# In this part we define the vector space and the torsion and curvature of the curve
print 
vs_base = UniformLagrangeVectorSpace(9)
vs = AffineVectorSpace(vs_base)
x = np.linspace(0,1,1025)

# these are parameters to write the curve
r = lambda s : 0.01 #s/(np.pi*8) + 1/np.pi*4
c = lambda t : t*0.1 + 0.1
T = 2 * np.pi

# s in the parameter of the curve and t the time
kappa =  lambda t : lambda s :  10*(1+s) #r(t)/(r(t)**2 + c(t)**2)
tau =  lambda t : lambda s : 10*(1.2+np.sin(2*T*t)) #c(t)/(r(t)**2 + c(t)**2)

curve = lambda t : ALCFromKappaAndTau( \
    vs, kappa(t), tau(t), s_space = np.linspace(0,1,2**7+1))

########
# PLOT #
########

xx = np.array([])
yy = np.array([])
zz = np.array([])

fig = plt.figure()
ax = Axes3D(fig)

ax.set_xlim3d([0.0, 1.0])
ax.set_xlabel('X')

ax.set_ylim3d([-.25, .25])
ax.set_ylabel('Y')

ax.set_zlim3d([-.25, .25])
ax.set_zlabel('Z')

line, = ax.plot(xx, yy, zz,  c="red", alpha=1)
head, = ax.plot([0], [0], [0],  marker='o', c="red", alpha=1)

# init is the starting point of the curve
def init():
    xx = np.array([])
    yy = np.array([])
    zz = np.array([])

# bar will print a progress bar in order to check the status of the output
bar = all_line_progress_bar()
# animate provide a frame for everey i in frames
def animate(i):
    bar.bar(i+1,nframes)
    xx=curve(float(i)/nframes).gamma(x)[0]
    yy=curve(float(i)/nframes).gamma(x)[1]
    zz=curve(float(i)/nframes).gamma(x)[2]
    ax.view_init(30, 45)
    line.set_data(xx,yy)
    line.set_3d_properties(zz)
    head.set_data(xx[0],yy[0])
    head.set_3d_properties(zz[0])

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=nframes, interval=20, blit=True)
anim.save('animation.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

print  "Done!\n"