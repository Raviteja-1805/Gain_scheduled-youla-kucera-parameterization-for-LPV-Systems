import sys
import numpy as np
import pickle
from design import design_lqr, design_lmi
from controller import check_bezout, ctrl_mats
from simulate import simulate, metrics, plot_results


method = sys.argv[1] if len(sys.argv) > 1 else "lqr"
sched = design_lmi() if method == "lmi" else design_lqr()

for th in (0.0, 0.5, 1.0):
    print(f"theta={th} Bezout errors:", check_bezout(sched, th))

Ak, Bk, Ck, H = ctrl_mats(0.5, sched, aw=True)
print("AW stable-factor poles @0.5:", np.linalg.eigvals(Ak - H @ Ck).round(3))

lin = simulate(sched, aw=False, umax=None)
no = simulate(sched, aw=False, umax=1.0)
yk = simulate(sched, aw=True, umax=1.0)


pickle.dump((lin, no, yk), open(f"logs_{method}.pkl", "wb"))


for name, lg in [("linear", lin), ("no AW", no), ("YK-AW", yk)]:
    print(name, metrics(lg))

plot_results(lin, no, yk, title=f"gains: {method}")