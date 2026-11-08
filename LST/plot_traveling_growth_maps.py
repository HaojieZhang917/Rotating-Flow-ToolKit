import csv
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent/'results'/(sys.argv[1] if len(sys.argv) > 1 else 'traveling_growth_maps')
f = lambda r, k: float(r[k])
best = list(csv.DictReader((root/'map_least_damped.csv').open()))
neut = list(csv.DictReader((root/'neutral_points.csv').open())) if (root/'neutral_points.csv').stat().st_size else []
# Thomas (2024) Fig. 2 styles
style = {0.0: ('tab:blue', '-'), -5.0: ('tab:red', '--'), 4.0: ('goldenrod', '-.'), 7.9: ('purple', ':'), 10.0: ('tab:green', '-')}
ws = [w for w in style if any(f(r, 'omega1') == w for r in best)]
IMAG_BB = 0.03   # points without a physical mode below this damping are shown as damped


def grid(w):
    s = [r for r in best if f(r, 'omega1') == w]
    allR = sorted({f(r, 'R') for r in best}); allm = sorted({f(r, 'm') for r in best})
    ms = allm
    G = np.full((len(ms), len(allR)), -IMAG_BB)            # -Im alpha; default: damped beyond the test window
    for r in s:
        G[ms.index(f(r, 'm')), allR.index(f(r, 'R'))] = -f(r, 'alpha_imag')
    return np.array(allR), np.array(ms), G


fig, axes = plt.subplots(2, 3, figsize=(16, 9), dpi=150)
for a, w in zip(axes.flat, ws):
    R, m, G = grid(w)
    c = a.pcolormesh(R, m, G, shading='nearest', cmap='RdBu_r', vmin=-0.1, vmax=0.1)
    a.contour(R, m, G, levels=[0], colors='k', linewidths=0.9)
    p = [r for r in neut if f(r, 'omega1') == w]
    a.plot([f(r, 'R') for r in p], [f(r, 'm') for r in p], '.', color='k', ms=2.5, label='refined neutral points')
    a.set_xlim(0, 800); a.set_ylim(-14, 80)
    a.set_title(rf'$\bar\omega=\omega R={w:g}$'); a.set_xlabel('R'); a.set_ylabel('n (= m)')
    fig.colorbar(c, ax=a, label=r'$-\mathrm{Im}\,\alpha$', extend='both')
    a.legend(frameon=False, fontsize=7, loc='upper left')
for a in list(axes.flat)[len(ws):]:
    a.axis('off')
fig.suptitle(r'Spatial growth rate $-\mathrm{Im}\,\alpha$ of the least-damped physical mode (rotating frame); black: neutral')
fig.tight_layout()
fig.savefig(root/'traveling_growth_maps.png'); fig.savefig(root/'traveling_growth_maps.pdf'); plt.close(fig)

fig, a = plt.subplots(figsize=(7.2, 5.4), dpi=170)
for w in ws:
    col, ls = style[w]
    R, m, G = grid(w)
    cs = a.contour(R, m, G, levels=[0], colors=[col], linestyles=[ls], linewidths=1.8)
    p = [r for r in neut if f(r, 'omega1') == w]
    a.plot([f(r, 'R') for r in p], [f(r, 'm') for r in p], '.', color=col, ms=2.5)
    a.plot([], [], color=col, ls=ls, lw=1.8, label=rf'$\bar\omega={w:g}$')
a.axhline(0, color='k', lw=0.6, ls=':')
a.set_xlim(0, 800); a.set_ylim(-12, 80); a.set_xlabel('Re (= R)'); a.set_ylabel('n')
a.set_title('Neutral curves in the (Re, n) plane (lines: grid zero contour; dots: refined)')
a.legend(frameon=True, fontsize=9, loc='upper left'); a.grid(True, alpha=.25)
fig.tight_layout()
fig.savefig(root/'neutral_curves_thomas_fig2.png'); fig.savefig(root/'neutral_curves_thomas_fig2.pdf'); plt.close(fig)
print('wrote', root)
