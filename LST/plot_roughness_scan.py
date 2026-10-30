import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out=os.path.join(os.path.dirname(__file__),"results","roughness_wavenumber_scan")
def read_csv(name):
    with open(os.path.join(out,name),newline="") as f:
        return list(csv.DictReader(f))
rows=read_csv("roughness_scan_main.csv")
k=np.array([float(r["kappa"]) for r in rows])
def arr(name):
    return np.array([float(r[name]) for r in rows])
note="R=285.36, alpha_f=0.38482, beta_f=0.07759, omega_f=-0.07759, N=109, zmax=30, z_ref=10, A_inc=1e-3, epsilon=1e-3"

def save(fig,name):
    fig.savefig(os.path.join(out,name+".png"),dpi=220)
    fig.savefig(os.path.join(out,name+".pdf"))
    plt.close(fig)

fig,ax=plt.subplots(figsize=(8,5))
for key,label,style in [("Gv_D","G_vol Dirichlet","--"),("Gw_D","G_wall Dirichlet","--"),("Gf_D","G_full Dirichlet","--"),("Gv_M","G_vol Modal","-"),("Gw_M","G_wall Modal","-"),("Gf_M","G_full Modal","-")]:
    ax.plot(k,arr(key),style,lw=1.8,label=label)
ax.set_yscale("log"); ax.set_xlabel(r"roughness radial wavenumber $\kappa$"); ax.set_ylabel(r"normalized gain $\mathcal{G}$")
ax.set_title("Second-order scattering gains versus roughness wavenumber")
ax.grid(True,which="both",alpha=.25); ax.legend(fontsize=8,ncol=2); fig.text(.02,.01,note,fontsize=7)
save(fig,"figure1_scattering_gains")

fig,axs=plt.subplots(3,1,figsize=(8,9),sharex=True)
axs[0].plot(k,arr("Rwall_D"),"--",label="D"); axs[0].plot(k,arr("Rwall_M"),"-",label="M"); axs[0].set_yscale("log"); axs[0].set_ylabel(r"$E_{wall}/E_{vol}$"); axs[0].legend(); axs[0].grid(alpha=.25)
axs[1].plot(k,arr("Ec_D")/arr("Ef_D"),"--",label="D"); axs[1].plot(k,arr("Ec_M")/arr("Ef_M"),"-",label="M"); axs[1].set_ylabel(r"$E_{cross}/E_{full}$"); axs[1].legend(); axs[1].grid(alpha=.25)
axs[2].plot(k,arr("rho_D"),"--",label="D"); axs[2].plot(k,arr("rho_M"),"-",label="M"); axs[2].set_xlabel(r"$\kappa$"); axs[2].set_ylabel(r"$\rho$"); axs[2].legend(); axs[2].grid(alpha=.25)
fig.suptitle("Relative wall contribution and interference",y=.995); fig.text(.02,.005,note,fontsize=7); fig.tight_layout(); save(fig,"figure2_coupling_contributions")

fig,ax=plt.subplots(figsize=(8,5))
for key,label in [("delta_vol",r"$\delta_{vol}$"),("delta_wall",r"$\delta_{wall}$"),("delta_full",r"$\delta_{full}$")]:
    ax.plot(k,arr(key),lw=1.8,label=label)
ax.set_yscale("log"); ax.set_xlabel(r"$\kappa$"); ax.set_ylabel("relative D/M difference on [0,10]"); ax.set_title("Dirichlet versus Modal second-order differences"); ax.grid(True,which="both",alpha=.25); ax.legend(); fig.text(.02,.01,note,fontsize=7); save(fig,"figure3_dirichlet_modal_difference")

fig,axs=plt.subplots(2,1,figsize=(8,8),sharex=True)
axs[0].plot(k,arr("source_D"),"--",label=r"$||S_{vol}||$ D"); axs[0].plot(k,arr("source_M"),"-",label=r"$||S_{vol}||$ M"); axs[0].set_yscale("log"); axs[0].set_ylabel("source norm"); axs[0].legend(); axs[0].grid(True,which="both",alpha=.25)
axs[1].plot(k,arr("Gv_D")/arr("source_D")**2,"--",label=r"$G_{vol}/||S||^2$ D"); axs[1].plot(k,arr("Gv_M")/arr("source_M")**2,"-",label=r"$G_{vol}/||S||^2$ M"); axs[1].set_yscale("log"); axs[1].set_xlabel(r"$\kappa$"); axs[1].set_ylabel("descriptive response/source ratio"); axs[1].legend(); axs[1].grid(True,which="both",alpha=.25)
fig.suptitle("Volume forcing strength and descriptive response amplification",y=.995); fig.text(.02,.005,note,fontsize=7); fig.tight_layout(); save(fig,"figure4_source_and_response")

prof=read_csv("roughness_scan_profiles.csv")
ks=sorted(set(float(r["kappa"]) for r in prof))
names=[("u","uD_re","uD_im","uM_re","uM_im"),("v","vD_re","vD_im","vM_re","vM_im"),("w","wD_re","wD_im","wM_re","wM_im")]
fig,axs=plt.subplots(1,3,figsize=(14,4),sharex=True)
for ax,(nm,dr,di,mr,mi) in zip(axs,names):
    for kval in ks:
        rr=[r for r in prof if float(r["kappa"])==kval]; z=np.array([float(r["z"]) for r in rr])
        d=np.sqrt(np.array([float(r[dr]) for r in rr])**2+np.array([float(r[di]) for r in rr])**2)
        m=np.sqrt(np.array([float(r[mr]) for r in rr])**2+np.array([float(r[mi]) for r in rr])**2)
        ax.plot(z,d,"--",label=f"D κ={kval:g}"); ax.plot(z,m,"-",label=f"M κ={kval:g}")
    ax.set_title(nm); ax.set_xlabel("z"); ax.grid(alpha=.25)
axs[0].set_ylabel(r"$|q_{hf}|$"); axs[-1].legend(fontsize=7,ncol=2); fig.suptitle("Typical second-order velocity profiles",y=1.02); fig.text(.02,-.02,note,fontsize=7); fig.tight_layout(); save(fig,"figure5_typical_profiles")

conv=read_csv("roughness_scan_convergence.csv")
fig,axs=plt.subplots(1,2,figsize=(11,4))
for kval in sorted(set(float(r["kappa"]) for r in conv)):
    rr=[r for r in conv if float(r["kappa"])==kval and float(r["zmax"])==30.0]
    N=np.array([float(r["N"]) for r in rr]); axs[0].plot(N,[float(r["Gf_D"]) for r in rr],"o--",label=f"D κ={kval:g}"); axs[0].plot(N,[float(r["Gf_M"]) for r in rr],"o-",label=f"M κ={kval:g}")
    rr=[r for r in conv if float(r["kappa"])==kval and float(r["N"])==109]
    L=np.array([float(r["zmax"]) for r in rr]); axs[1].plot(L,[float(r["Gf_D"]) for r in rr],"o--",label=f"D κ={kval:g}"); axs[1].plot(L,[float(r["Gf_M"]) for r in rr],"o-",label=f"M κ={kval:g}")
axs[0].set_yscale("log"); axs[1].set_yscale("log"); axs[0].set_xlabel("N at zmax=30"); axs[1].set_xlabel("zmax at N=109"); axs[0].set_ylabel("G_full"); axs[1].set_ylabel("G_full"); axs[0].grid(True,which="both",alpha=.25); axs[1].grid(True,which="both",alpha=.25); axs[1].legend(fontsize=7)
fig.suptitle("Representative spectral and outer-boundary convergence",y=1.02); fig.text(.02,-.02,note,fontsize=7); fig.tight_layout(); save(fig,"figure6_convergence")
print("plots written")
