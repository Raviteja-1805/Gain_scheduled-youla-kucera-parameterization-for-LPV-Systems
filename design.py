import numpy as np
import scipy.linalg as la
from plant import A_of, B, C, n, m, p, augment

try:

    import cvxpy as cp
    HAVE_CVXPY = True
except ImportError:
    HAVE_CVXPY = False


def lqr(Aa, Ba, Q, R):
    P = la.solve_continuous_are(Aa, Ba, Q, R)
    return -np.linalg.solve(R, Ba.T @ P)

def kalman(A, C, W, V):
    P = la.solve_continuous_are(A.T, C.T, W, V)
    return P @ C.T @ np.linalg.inv(V)

class GainSchedule:
    def __init__(self, grid, F_tab, L_tab):
        self.grid = np.asarray(grid, float)
        self.F_tab = np.asarray(F_tab)
        self.L_tab = np.asarray(L_tab)

    def _interp(self, th, tab):
        g = self.grid
        th = np.clip(th, g[0], g[-1])
        i = min(np.searchsorted(g, th, side='right') - 1, len(g) - 2)
        w = (th - g[i]) / (g[i + 1] - g[i])
        return (1- w) * tab[i] + w * tab[i + 1]

    def F(self, th): return self._interp(th, self.F_tab)
    def L(self, th): return self._interp(th, self.L_tab)

def design_lqr(n_grid=5):
    Qw, Rw = np.diag([4.0, 1.0, 6.0]), np.eye(1)
    W, V = np.eye(2), 0.1 * np.eye(1)
    grid = np.linspace(0, 1, n_grid)
    F_tab = [lqr(*augment(A_of(t), B, C), Qw, Rw) for t in grid]
    L_tab = [kalman(A_of(t), C, W, V) for t in grid]
    return GainSchedule(grid, F_tab, L_tab)

def design_lmi(alpha_F=0.6, alpha_L=2.0, eps=1e-4):
    if not HAVE_CVXPY:
        raise ImportError("design_lmi needs cvxpy")
    verts = [0.0, 1.0]
    na = n + p

    # ---------- state feedback on the augmented plant (3x3)
    P  = cp.Variable((na, na), symmetric=True)
    Ys = [cp.Variable((m, na)) for _ in verts]
    cons = [P >> np.eye(na)]
    for th, Y in zip(verts, Ys):
        Aa, Ba = augment(A_of(th), B, C)
        M = Aa @ P + P @ Aa.T + Ba @ Y + Y.T @ Ba.T + 2 * alpha_F * P
        cons.append((M + M.T) / 2 << -eps * np.eye(na))
    cp.Problem(cp.Minimize(sum(cp.norm(Y, 'fro') for Y in Ys)), cons).solve()
    F_v = [Y.value @ np.linalg.inv(P.value) for Y in Ys]

    # ---------- observer (2x2)
    Q  = cp.Variable((n, n), symmetric=True)
    Ws = [cp.Variable((n, p)) for _ in verts]
    cons = [Q >> np.eye(n)]
    for th, Wv in zip(verts, Ws):
        A = A_of(th)
        M = A.T @ Q + Q @ A - Wv @ C - C.T @ Wv.T + 2 * alpha_L * Q
        cons.append((M + M.T) / 2 << -eps * np.eye(n))
    cp.Problem(cp.Minimize(sum(cp.norm(Wv, 'fro') for Wv in Ws)), cons).solve()
    L_v = [np.linalg.solve(Q.value, Wv.value) for Wv in Ws]

    return GainSchedule(verts, F_v, L_v)