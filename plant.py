import numpy as np

A0 = np.array([[0., 1.],[-1, -0.4]])
A1 = np.array([[0., 0.], [0.9, 0.]])
B  = np.array([[0.], [1.]])
C  = np.array([[1., 0.]])
n, m, p = 2, 1, 1

def A_of(th):
    return A0 + th * A1

def augment(A, B, C):
    n, m, p = A.shape[0], B.shape[1], C.shape[0]
    Aa = np.block ([[A, np.zeros((n, p))], 
                    [-C, np.zeros((p, p))]])
    Ba = np.vstack([B, np.zeros((p, m))])
    return Aa, Ba

def th_fun(t):
    return 0.5 + 0.5 * np.sin(2 * np.pi * t/ 40)

def r_fun(t):
    return 2.0 if 2 <= t < 22 else (-1.5 if 22 <= t < 42 else 0.5)