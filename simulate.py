import numpy as np
import matplotlib.pyplot as plt
from plant import A_of, B, C, n, m, p, th_fun, r_fun
from controller import ctrl_mats

def simulate(sched, aw, umax, T=60.0, dt=1e-3, kaw=6.0):
    N = int(T / dt)
    z = np.zeros(n + n + p)
    log = {k: np.zeros(N) for k in ('t', 'y', 'u', 'xi', 'th', 'r')}

    def sat(v):
        return np.clip(v, -umax, umax) if umax is not None else v

    def f(t, z):
        x, xk = z[:n], z[n:]
        th, r = th_fun(t), np.atleast_1d(r_fun(t))
        Ak, Bk, Ck, H = ctrl_mats(th, sched, aw, kaw)
        y, uc = C @ x, Ck @ xk
        u = sat(uc)
        dx = A_of(th) @ x + B @ u
        dxk = Ak @ xk + Bk @ np.concatenate([y, r]) + H @ (u - uc)
        return np.concatenate([dx, dxk])

    for k in range(N):
        t = k * dt
        x, xk = z[:n], z[n:]
        th = th_fun(t)
        uc = ctrl_mats(th, sched, aw, kaw)[2] @ xk
        log['t'][k], log['y'][k], log['u'][k] = t, (C @ x)[0], sat(uc)[0]
        log['xi'][k], log['th'][k], log['r'][k] = xk[n], th, r_fun(t)
        k1 = f(t, z);
        k2 = f(t + dt/2, z + dt/2 * k1)
        k3 = f(t + dt/2, z + dt/2 * k2)
        k4 = f(t + dt, z + dt * k3)
        z = z + dt / 6 * (k1 + 2*k2 + 2*k3 + k4 )
    return log

def metrics(lg):
    dt = lg['t'][1] - lg['t'][0]
    return np.sum(abs(lg['r'] - lg['y'])) * dt, np.max(abs(lg['xi']))

def plot_results(lin, no, yk, title="", method="lqr"):
    tag = method.upper()
    fig, ax = plt.subplots(3, 1, figsize=(9, 8), sharex=True)

    ax[0].plot(yk['t'], yk['r'], 'k--', label='Reference')
    ax[0].plot(lin['t'], lin['y'], 'g:', label='Linear (no saturation)')
    ax[0].plot(no['t'],  no['y'],  'r',  label='Windup (no anti-windup)')
    ax[0].plot(yk['t'],  yk['y'],  'b',  label=f'YK anti-windup ({tag})')
    ax[0].set_ylabel('y'); ax[0].set_title(title)

    ax[1].plot(no['t'], no['u'], 'r', label='Windup')
    ax[1].plot(yk['t'], yk['u'], 'b', label=f'YK-AW ({tag})')
    ax[1].set_ylabel('u')

    ax[2].plot(no['t'], no['xi'], 'r', label='Windup')
    ax[2].plot(yk['t'], yk['xi'], 'b', label=f'YK-AW ({tag})')
    ax[2].set_ylabel('integrator'); ax[2].set_xlabel('t [s]')

    for a in ax:
        a.grid(alpha=.3)
        a.legend(loc='upper right', fontsize=7, frameon=True,
                 framealpha=0.9, borderpad=0.4, handlelength=1.5)

    plt.tight_layout()
    plt.savefig(f"results_{method}.png", dpi=130)
    print(f"saved results_{method}.png")
    plt.show()