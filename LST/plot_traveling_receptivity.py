import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent/'results'/'traveling_receptivity'
rows = list(csv.DictReader((root/'traveling_receptivity_scan.csv').open()))
f = lambda r, k: float(r[k])
titles = {'TypeII_w8': r'Type-II, $\bar\omega=8$, $R=200$, $m_t=-2$ ($m_f=-8$, $n=6$)',
          'TypeI_w4': r'slow Type-I, $\bar\omega=4$, $R=400$, $m_t=22$ ($m_f=-4$, $n=26$)'}
colors = {'poloidal': 'tab:blue', 'toroidal': 'tab:orange'}
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), dpi=160)
for a, (tg, title) in zip(ax, titles.items()):
    for pol, c in colors.items():
        rs = sorted([r for r in rows if r['target'] == tg and r['polarization'] == pol], key=lambda r: f(r, 'alpha_f'))
        for sign in (-1, 1):
            s = [r for r in rs if np.sign(f(r, 'alpha_f')) == sign]
            af = np.array([f(r, 'alpha_f') for r in s]); L = np.array([f(r, 'abs_Lambda') for r in s])
            ok = np.array([f(r, 'rho_farfield') < 1 for r in s])
            a.semilogy(af, np.where(ok, L, np.nan), 'o-', color=c, ms=3, lw=1.3, label=pol if sign > 0 else None)
            a.semilogy(af[~ok], L[~ok], 'x', color=c, ms=5, alpha=0.6,
                       label=(r'$\rho\geq1$ (upper-deck divergent)' if (sign > 0 and pol == 'poloidal') else None))
        qh = np.array([f(r, 'qh_norm') for r in rs]); afa = np.array([f(r, 'alpha_f') for r in rs])
    for x0 in afa[qh > 4]:
        a.axvspan(x0 - 0.025, x0 + 0.025, color='tab:red', alpha=0.08, lw=0)
    a.set_title(title, fontsize=10); a.set_xlabel(r'$\alpha_f$  ($\kappa=\alpha_n-\alpha_f$)')
    a.set_ylabel(r'$|\Lambda|$ (edge-normalised vortex, unit roughness, $v(z^*)=1$)')
    a.grid(True, which='both', alpha=.3); a.legend(frameon=False, fontsize=8)
fig.suptitle(r'Vortex–roughness receptivity of traveling modes (complex-$\alpha$ residue); red band: $\|q_h\|>4$ (stationary self-resonance of the roughness)')
fig.tight_layout()
fig.savefig(root/'traveling_receptivity_scan.png'); fig.savefig(root/'traveling_receptivity_scan.pdf'); plt.close(fig)
print('wrote', root)
