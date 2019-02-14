import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

sig1 = 1
sig2 = 2
theta = np.pi/4

SigmaU =  np.array([[sig1, 0], [0, sig2]])

R = np.array([[np.cos(theta), -np.sin(theta)],
              [np.sin(theta), np.cos(theta)]])

SigmaX = R@SigmaU@R.T

ax = plt.subplot(121)

vals, vecs = eigsorted(SigmaU)
th = np.degrees(np.arctan2(*vecs[:,0][::-1]))
w, h = np.sqrt(vals)
for i in range(1,4):
    ell = Ellipse(xy=(0, 0),
                  width=i*w, height=i*h,
                  angle=th, color='black')
    ell.set_facecolor('none')
    ax.add_artist(ell)
ax.set_xlim([-3,3])
ax.set_ylim([-3,3])
ax.set_title("$U_1, U_2$ curves of constant probability")
ax.set_xlabel("$u_1$")
ax.set_ylabel("$u_2$")

ax = plt.subplot(122)
vals, vecs = eigsorted(SigmaX)
th = np.degrees(np.arctan2(*vecs[:,0][::-1]))
w, h = np.sqrt(vals)
for i in range(1,4):
    ell = Ellipse(xy=(0, 0),
                  width=i*w, height=i*h,
                  angle=th, color='black')
    ell.set_facecolor('none')
    ax.add_artist(ell)
ax.set_xlim([-3,3])
ax.set_ylim([-3,3])
ax.set_title("$X_1, X_2$ curves of constant probability")
ax.set_xlabel("$x_1$")
ax.set_ylabel("$x_2$")

plt.show()
