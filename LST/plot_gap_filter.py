"""Figure for the gap-filter verification (analyze_gap_filter*.jl, export_gap_filter_profiles.jl)."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent / 'results' / 'gap_filter'
f = lambda p: [{k: (float(v) if k != 'case' else v) for k, v in r.items()} for r in csv.DictReader((root / p).open())]
core, E, prof = f('core_test.csv'), f('partE.csv'), f('profiles.csv')
INK, MUTED = '#2b2b2b', '#8a8a85'
plt.rcParams.update({'axes.grid': True, 'grid.color': '#e4e3de', 'font.size': 9})
fig, ax = plt.subplots(1, 4, figsize=(17, 4.4), dpi=150)

# (a) core decay: numerical vs theory
ev = [r for r in core if abs(r['sig_over_f']) > 2]; it = [r for r in core if abs(r['sig_over_f']) < 0.8]
x = [r['lnratio_inviscid'] for r in ev]; y = [r['lnratio_num'] for r in ev]
ax[0].plot(x, y, 'o', ms=3.5, color='#2a78d6', label=r'evanescent $|\sigma/f|>2$ (240 cases)')
ax[0].plot([r['lnratio_inviscid'] for r in it], [r['lnratio_num'] for r in it], 's', ms=3.5, mfc='none', color='#eb6834',
           label=r'inertial $|\sigma/f|<0.8$ (64, all with critical level)')
lim = [min(x + y) - 0.5, max(x + y) + 0.5]; ax[0].plot(lim, lim, color=INK, lw=1)
ax[0].set_xlabel(r'theory  $\ln|\cosh m(\zeta_b-\zeta_r)/\cosh m(\zeta_a-\zeta_r)|$'); ax[0].set_ylabel(r'numerical  $\ln|u(z_b)/u(z_a)|$')
ax[0].set_title('(a) core decay, z = 0.2 → 0.5 (inviscid m, rigid lid)'); ax[0].legend(fontsize=7, frameon=False)

# (b) direct path vs Re_h
md = np.array([r['m'] * r['Dzeta'] for r in E]); lL = np.array([r['lnLdir'] for r in E])
ax[1].plot(md, lL, 'o', color='#2a78d6', ms=6, label=r'$\ln|\Lambda_{dir}|$, $Re_h$=1000…4000')
xx = np.linspace(md.min() - 1, md.max() + 1, 10); ax[1].plot(xx, -xx - 0.48, color=INK, lw=1, label=r'$-m\Delta\zeta-0.48$ (slope −1)')
for r, xv, yv in zip(E, md, lL):
    ax[1].annotate(f"{int(r['Re_h'])}", (xv, yv), textcoords='offset points', xytext=(5, 3), fontsize=7, color=MUTED)
ax[1].set_xlabel(r'$m\,\Delta\zeta$'); ax[1].set_ylabel(r'$\ln|\Lambda_{dir}|$'); ax[1].legend(fontsize=7, frameon=False)
ax[1].set_title('(b) direct path, Type-I R=330, n=16, ω̄=8')

# (c)/(d) profiles
for k, (case, title) in enumerate((('evanescent_a0.45_n16_wbar8', r'(c) evanescent: $\alpha$=0.45, n=16, ω̄=8 ($\sigma/f$=30)'),
                                     ('inertial_a0.2_n-10_c0.7', r'(d) inertial: $\alpha$=0.2, n=−10, c=0.7 ($\sigma_{core}/f$=0.23)'))):
    a = ax[2 + k]; q = [r for r in prof if r['case'] == case]
    z = np.array([r['z'] for r in q]); u = np.array([r['absu'] for r in q]); sg = np.array([r['sigma_local_over_f'] for r in q])
    a.semilogx(u / u[np.argmin(abs(z - 0.5))], z, color='#2a78d6', lw=1.6, label='|u| (numerical, normalised at z=0.5)')
    if k == 0:
        zz = np.linspace(0.124, 0.85, 50); m = 0.4529; s = np.sqrt(1000)
        a.semilogx(np.cosh(m * (zz - 0.124) * s) / np.cosh(m * (0.5 - 0.124) * s), zz, '--', color=INK, lw=1.2, label=r'theory $\cosh m(\zeta-\zeta_r)$')
    else:
        zc = [0.5 * (z[i] + z[i + 1]) for i in range(len(z) - 1) if 0.12 < z[i] < 0.9 and np.sign(sg[i]) != np.sign(sg[i + 1])]
        for j, zc_ in enumerate(zc):
            a.axhline(zc_, color='#eb6834', lw=0.9, ls=':', label='critical level σ(z)=0' if j == 0 else None)
    for yv in (0.124, 0.845):
        a.axhline(yv, color=MUTED, lw=0.7)
    a.set_xlabel('|u| / |u(0.5)|'); a.set_ylabel('z (rotor 0, stator 1)'); a.set_title(title, fontsize=9); a.legend(fontsize=7, frameon=False, loc='upper left')
fig.suptitle('Gap filter: verification of the rotating-core transfer theory ($Re_h$=1000 unless stated; grey lines: layer edges z_r, z_s)')
fig.tight_layout(); fig.savefig(root / 'gap_filter_verification.png'); fig.savefig(root / 'gap_filter_verification.pdf'); plt.close(fig)
print('wrote', root)
