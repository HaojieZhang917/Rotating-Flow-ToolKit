from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
TAG = r'$R=285.36$, $(\\alpha_f,\\beta_f,\\omega_f)=(0.38482,0.07759,-0.07759)$, $A_{inc}=\\epsilon=10^{-3}$'

def read_csv(name):
    with (ROOT / name).open() as f:
        return list(csv.DictReader(f))

def save(fig, name):
    fig.tight_layout()
    fig.savefig(ROOT / (name + '.png'), dpi=180)
    fig.savefig(ROOT / (name + '.pdf'))
    plt.close(fig)

ft = read_csv('fourier_transform_validation.csv')
sp = read_csv('roughness_spectrum.csv')
r = np.array([float(x['r']) for x in ft])
he = np.array([float(x['h_exact']) for x in ft])
hr = np.array([float(x['h_reconstructed']) for x in ft])
k = np.array([float(x['kappa']) for x in sp])
hh = np.array([float(x['abs_Hhat']) for x in sp])

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(r, he, label='analytic Gaussian h/e')
ax[0].plot(r, hr, '--', label='inverse Fourier reconstruction')
ax[0].set(xlabel=r'$r-r_f$', ylabel=r'$h(r)/\\epsilon$', title='Gaussian roughness geometry')
ax[0].legend(frameon=False); ax[0].grid(alpha=.3)
ax[1].plot(k, hh)
ax[1].set(xlabel=r'$\\kappa$ (radial Fourier wavenumber)', ylabel=r'$|\\hat H(\\kappa)|$', title='Gaussian Fourier spectrum')
ax[1].grid(alpha=.3)
fig.suptitle('Figure 1: roughness model and Fourier spectrum\\n' + TAG)
save(fig, 'figure1_roughness_model')

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.semilogy(r, np.maximum(abs(he - hr), 1e-16))
ax.set(xlabel=r'$r-r_f$', ylabel='absolute inverse-transform error', title='Figure 2: Fourier reconstruction error')
ax.grid(alpha=.3); ax.text(.02, .03, TAG, transform=ax.transAxes, fontsize=8)
save(fig, 'figure2_fourier_validation')

g = np.array([float(x['G_single']) for x in sp])
wg = np.array([float(x['abs_Hhat_sq_G_single']) for x in sp])
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(k, g, 'o-', ms=3, label='unit Fourier component')
ax.axvline(-.185, color='tab:red', ls='--', label=r'$\\kappa=-0.185$')
ax.axvline(.185, color='tab:green', ls='--')
ax.set(xlabel=r'$\\kappa$ (signed roughness wavenumber)', ylabel=r'$G_{single}=E_{hf}/|A_{inc}|^2$', title='Figure 3: single-Fourier scattering transfer')
ax.legend(frameon=False); ax.grid(alpha=.3); save(fig, 'figure3_single_mode_gain')

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(k, wg, 'o-', ms=3, color='tab:purple')
ax.set(xlabel=r'$\\kappa$', ylabel=r'$|\\hat H(\\kappa)|^2G_{single}(\\kappa)$', title='Figure 4: Gaussian-weighted scattering spectrum')
ax.grid(alpha=.3); ax.text(.02, .03, 'Spectral weight only; not local energy at fixed r', transform=ax.transAxes, fontsize=8)
save(fig, 'figure4_weighted_spectrum')

rows = read_csv('reconstructed_response.csv')
rv = np.array([float(x['r']) for x in rows])
zv = np.array([float(x['z']) for x in rows])
rs = np.unique(rv); zs = np.unique(zv)
shape = (len(rs), len(zs))
def field(key):
    vals = np.array([complex(float(x[key + '_re']), float(x[key + '_im'])) for x in rows])
    return vals.reshape(shape)
u, v, w = field('u'), field('v'), field('w')
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for a, F, name in zip(ax, [u, v, w], ['$|u_{hf}|$', '$|v_{hf}|$', '$|w_{hf}|$']):
    im = a.pcolormesh(rs, zs, np.abs(F).T, shading='auto')
    a.set(xlabel=r'$r$ (local radial coordinate)', ylabel=r'$z$ (wall-normal coordinate)', title=name)
    fig.colorbar(im, ax=a)
fig.suptitle('Figure 5: complete reconstructed second-order response\\n' + TAG)
save(fig, 'figure5_reconstructed_response')

le = read_csv('local_energy.csv')
rr = np.array([float(x['r']) for x in le])
ee = np.array([float(x['E_hf']) for x in le])
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.plot(rr, ee)
ax.set(xlabel=r'$r$ (local radial coordinate)', ylabel=r'$E_{hf}(r)=\\int_0^{z_i}|\\boldsymbol u_{hf}|^2dz$', title='Figure 6: local reconstructed scattering energy')
ax.grid(alpha=.3); ax.text(.02, .03, 'Complex Fourier components are summed before forming energy', transform=ax.transAxes, fontsize=8)
save(fig, 'figure6_local_energy')

cv = read_csv('reconstruction_convergence.csv')
labels = [f"K={x['kappa_max']}, Nκ={x['nk']}" for x in cv]
ep = np.array([float(x['relative_to_main']) for x in cv])
geom = np.array([float(x['geometry_error']) for x in cv])
fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(len(labels))
ax.semilogy(x, np.maximum(ep, 1e-16), 'o-', label='response energy difference')
ax.semilogy(x, np.maximum(geom, 1e-16), 's--', label='Gaussian geometry error')
ax.set(xticks=x, xticklabels=labels, xlabel='spectral truncation and grid', ylabel='relative error', title='Figure 7: Fourier reconstruction convergence')
ax.tick_params(axis='x', rotation=45); ax.legend(frameon=False); ax.grid(alpha=.3)
save(fig, 'figure7_reconstruction_convergence')
