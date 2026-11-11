import numpy as np
from plant import A_of, B, C, n, m, p

def tf(A_, B_, C_, D_, s):
    return C_ @ np.linalg.solve(s * np.eye(A_.shape[0]) - A_, B_) + D_

def check_bezout(sched, th, w=0.5):
    A = A_of(th)
    Fx = sched.F(th)[:, :n]
    L = sched.L(th)
    s = 1j * w
    I, Z = np.eye(m), np.zeros((p, m))
    M = tf(A + B@Fx, B, Fx, I, s);
    N = tf(A + B@Fx, B, C, Z, s);
    Mt = tf(A - L@C, -L, C, np.eye(p), s);
    Nt = tf(A - L@C, B, C, Z, s);
    X = tf(A - L@C, -B, Fx, I, s);
    Y = tf(A - L@C, L, Fx, np.zeros((m, p)), s)
    G = tf(A, B, C, Z, s)
    return (abs(N @ np.linalg.inv(M) - G).max(),
            abs(np.linalg.inv(Mt) @ Nt - G).max(),
            abs(X @ M - Y @ N - I).max())

def ctrl_mats(th, sched, aw=True, kaw=6.0):
    A = A_of(th)
    Fa = sched.F(th)
    L = sched.L(th)
    Fx, Fi = Fa[:, :n], Fa[:, n:]
    
    

    Ak = np.block([[A - L@C + B@Fx, B@Fi],
                   [np.zeros((p, n)), np.zeros((p, p))]])
    Bk = np.block([[L, np.zeros((n, p))],
                   [-np.eye(p), np.eye(p)]])
    Ck = Fa
    if aw:
        H = np.vstack([B, kaw * np.linalg.inv(Fi)])
    else:
        H = np.zeros((n + p, m))
    return Ak, Bk, Ck, H