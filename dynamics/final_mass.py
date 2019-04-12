import numpy as np


def parallel_axis(I_G, m, r_B_G):
    x = r_B_G[0]
    y = r_B_G[1]
    z = r_B_G[2]

    I_add = np.array([[y**2 + z**2, -x*y, -x*z],
                      [-x*y, x**2 + z**2, -y*z],
                      [-x*z, -y*z, x**2 + y**2]])
    return I_G + m*I_add

def calc_com(masses, vectors):
    M = np.sum(masses)
    return np.sum(vectors*masses,axis=1)/M

def calc_inertia(com, masses, vectors, inertias):
    r_B_G = -vectors + com[:,np.newaxis]
    I = [parallel_axis(inertias[i], m, r_B_G[:,i]) for i,m in enumerate(masses)]
    return reduce(np.ndarray.__add__, I)

# Parameters
I2M = 1/39.3701
LB2KG = 0.45359237
SLUGFT2KGM = 1.355817


CM_M = 9730*LB2KG
SM_M = 9690*LB2KG
PROP_M = 37295*LB2KG

masses = np.array([CM_M, SM_M, PROP_M])

CM_I = np.diag([4474., 3919., 3684.])*SLUGFT2KGM
SM_I = np.diag([6222., 10321., 10136.])*SLUGFT2KGM
PROP_I = np.diag([19162., 19872., 26398.])*SLUGFT2KGM
inertias = [CM_I, SM_I, PROP_I]

# General Case
r_CM_A =  np.array([1043.1, -0.1, 7.8])*I2M
r_SM_A = np.array([908.2, 0.7, -0.6])*I2M
r_PROP_A = np.array([905.9, 5.6, -2.4])*I2M
vectors = np.column_stack([r_CM_A, r_SM_A, r_PROP_A])

r_CG_A = calc_com(masses, vectors)
CG_I_gen = calc_inertia(r_CG_A, masses, vectors, inertias)

# Specific Case
r_CM_A =  np.array([1043.1, 0., 0.])*I2M
r_SM_A = np.array([908.2, 0., 0.])*I2M
r_PROP_A = np.array([905.9, 0., 0.])*I2M
vectors = np.column_stack([r_CM_A, r_SM_A, r_PROP_A])

r_CG_A = calc_com(masses, vectors)
CG_I_sim = calc_inertia(r_CG_A, masses, vectors, inertias)
