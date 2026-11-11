import sys
import pickle
from simulate import plot_results

method = sys.argv[1] if len (sys.argv) > 1 else "lqr"
lin, no, yk = pickle.load(open(f"logs_{method}.pkl", "rb"))
plot_results(lin, no, yk, title=f"gains: {method}")