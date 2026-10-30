import csv
from pathlib import Path
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'second_order_modal_peak_scan.csv').open()))
k=[float(r['kappa']) for r in rows]
gd=[float(r['G_D']) for r in rows]
gm=[float(r['G_M']) for r in rows]
rel=[float(r['D_M_relative'])*100 for r in rows]

plt.figure(figsize=(6.4,4.5),dpi=180)
plt.plot(k,gd,'o-',lw=1.5,ms=3,label='six-state Dirichlet')
plt.plot(k,gm,'s-',lw=1.5,ms=3,label='six-state Modal')
plt.xlabel(r'$|\kappa|$ (signed roughness wave number $\kappa_s=-|\kappa|$)')
plt.ylabel(r'$G_{\mathrm{full}}=\|u_{hf}\|_{L^2(0,z_i)}^2/(\epsilon^2 |A_{inc}|^2)$')
plt.title('Second-order roughness scattering gain')
plt.grid(True,alpha=.3); plt.legend(frameon=False); plt.tight_layout()
plt.savefig(root/'second_order_modal_gain_scan.png'); plt.savefig(root/'second_order_modal_gain_scan.pdf'); plt.close()

plt.figure(figsize=(6.4,4.2),dpi=180)
plt.plot(k,rel,'o-',lw=1.5,ms=3,color='tab:red')
plt.xlabel(r'$|\kappa|$ (negative scattering branch)')
plt.ylabel(r'$100|G_M-G_D|/G_M$ (%)')
plt.title('Modal versus Dirichlet boundary-model difference')
plt.grid(True,alpha=.3); plt.tight_layout()
plt.savefig(root/'second_order_modal_relative_difference.png'); plt.savefig(root/'second_order_modal_relative_difference.pdf'); plt.close()
