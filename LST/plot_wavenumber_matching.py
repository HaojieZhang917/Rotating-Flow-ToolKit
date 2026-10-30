import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out=os.path.join(os.path.dirname(__file__),"results","wavenumber_matching")
def read(name):
    with open(os.path.join(out,name),newline="") as f: return list(csv.DictReader(f))
def a(rows,key): return np.array([float(r[key]) for r in rows])
def save(fig,name):
    fig.savefig(os.path.join(out,name+".png"),dpi=220); fig.savefig(os.path.join(out,name+".pdf")); plt.close(fig)
note="beta_f=0.07759, omega_f=-0.07759, A_inc=1e-3, epsilon=1e-3, energy on [0,10], second-order outer BC=homogeneous Dirichlet"

low=read("small_kappa_scan.csv")
fig,axs=plt.subplots(3,1,figsize=(8,9),sharex=True)
k=a(low,"kappa")
axs[0].plot(k,a(low,"Gf"),"o-",label="G_full"); axs[0].plot(k,a(low,"Gv"),"o--",label="G_vol"); axs[0].plot(k,a(low,"Gw"),"o--",label="G_wall"); axs[0].set_yscale("log"); axs[0].set_ylabel("gain"); axs[0].legend(); axs[0].grid(True,which="both",alpha=.25)
axs[1].plot(k,a(low,"source_norm"),"o-",label="||S_vol||"); axs[1].set_yscale("log"); axs[1].set_ylabel("source norm"); axs[1].legend(); axs[1].grid(True,which="both",alpha=.25)
axs[2].plot(k,a(low,"H_eff"),"o-",label="H_eff"); axs[2].set_xlabel(r"$\kappa$"); axs[2].set_ylabel(r"$\sqrt{E_{vol}}/||S_{vol}||$"); axs[2].legend(); axs[2].grid(True,alpha=.25)
fig.suptitle("Low-kappa limit",y=.995); fig.text(.02,.005,note,fontsize=7); fig.tight_layout(); save(fig,"figure1_low_kappa")

br=read("positive_negative_branches.csv")
fig,axs=plt.subplots(2,1,figsize=(8,7),sharex=True)
for b,sty in [("plus","-"),("minus","--")]:
    rr=[r for r in br if r["branch"]==b]; rr=sorted(rr,key=lambda r:float(r["kappa"]))
    axs[0].plot(a(rr,"kappa"),a(rr,"Gf"),sty,marker="o",label=b); axs[1].plot(a(rr,"kappa"),a(rr,"Gv"),sty,marker="o",label=b)
axs[0].set_yscale("log"); axs[1].set_yscale("log"); axs[0].set_ylabel("G_full"); axs[1].set_ylabel("G_vol"); axs[1].set_xlabel(r"positive physical $|\kappa|$"); axs[0].legend(); axs[1].legend(); axs[0].grid(True,which="both",alpha=.25); axs[1].grid(True,which="both",alpha=.25)
fig.suptitle("Positive and negative Fourier scattering branches",y=.995); fig.text(.02,.005,note,fontsize=7); fig.tight_layout(); save(fig,"figure2_positive_negative")

ar=read("alpha_kappa_scan.csv")
alphas=sorted(set(float(r["alpha_f"]) for r in ar)); kappas=sorted(set(float(r["kappa"]) for r in ar)); A=np.full((len(alphas),len(kappas)),np.nan); S=A.copy(); H=A.copy()
for r in ar:
    i=alphas.index(float(r["alpha_f"])); j=kappas.index(float(r["kappa"])); A[i,j]=float(r["Gf"]); S[i,j]=float(r["source_norm"]); H[i,j]=float(r["H_eff"])
fig,ax=plt.subplots(figsize=(8,5)); im=ax.pcolormesh(kappas,alphas,np.log10(A),shading="auto",cmap="viridis"); fig.colorbar(im,ax=ax,label=r"$\log_{10} G_{full}$"); ax.set_xlabel(r"$\kappa$"); ax.set_ylabel(r"$\alpha_f$"); ax.set_title("Two-dimensional positive-branch scattering gain"); fig.text(.02,.01,note,fontsize=7); fig.tight_layout(); save(fig,"figure3_alpha_kappa_gain")

fig,axs=plt.subplots(1,2,figsize=(12,4))
im=axs[0].pcolormesh(kappas,alphas,np.log10(S),shading="auto",cmap="magma"); fig.colorbar(im,ax=axs[0],label=r"$\log_{10}||S_{vol}||$"); axs[0].set_xlabel(r"$\kappa$"); axs[0].set_ylabel(r"$\alpha_f$")
im=axs[1].pcolormesh(kappas,alphas,np.log10(H),shading="auto",cmap="plasma"); fig.colorbar(im,ax=axs[1],label=r"$\log_{10} H_{eff}$"); axs[1].set_xlabel(r"$\kappa$"); axs[1].set_ylabel(r"$\alpha_f$")
fig.suptitle("Source strength and descriptive effective response",y=1.02); fig.text(.02,-.02,note,fontsize=7); fig.tight_layout(); save(fig,"figure4_source_heff")

peaks=read("coarse_peaks.csv"); ap=a(peaks,"alpha_f"); kp=a(peaks,"kappa"); at=a(peaks,"target_alpha")
fig,ax=plt.subplots(figsize=(8,5)); ax.plot(ap,kp,"o-",label=r"$\kappa_{max}$"); ax2=ax.twinx(); ax2.plot(ap,at,"s--",color="tab:red",label=r"$\alpha_{hf,max}$"); ax.set_xlabel(r"$\alpha_f$"); ax.set_ylabel(r"$\kappa_{max}$"); ax2.set_ylabel(r"$\alpha_{hf,max}$",color="tab:red"); ax.set_title("Location of coarse-scan maxima"); ax.grid(alpha=.25); fig.text(.02,.01,note,fontsize=7); fig.tight_layout(); save(fig,"figure5_peak_matching")

prof_path=os.path.join(os.path.dirname(__file__),"results","roughness_wavenumber_scan","roughness_scan_profiles.csv")
if os.path.exists(prof_path):
    with open(prof_path,newline="") as f: pr=list(csv.DictReader(f))
    ks=sorted(set(float(r["kappa"]) for r in pr)); names=[("u","uD_re","uD_im","uM_re","uM_im"),("v","vD_re","vD_im","vM_re","vM_im"),("w","wD_re","wD_im","wM_re","wM_im")]
    fig,axs=plt.subplots(1,3,figsize=(14,4),sharex=True)
    for ax,(nm,dr,di,mr,mi) in zip(axs,names):
        for kval in ks:
            rr=[r for r in pr if float(r["kappa"])==kval]; z=np.array([float(r["z"]) for r in rr]); d=np.hypot(a(rr,dr),a(rr,di)); m=np.hypot(a(rr,mr),a(rr,mi)); ax.plot(z,d,"--",label=f"D {kval:g}"); ax.plot(z,m,"-",label=f"M {kval:g}")
        ax.set_title(nm); ax.set_xlabel("z"); ax.grid(alpha=.25)
    axs[0].set_ylabel(r"$|q_{hf}|$"); axs[-1].legend(fontsize=7,ncol=2); fig.suptitle("Representative response profiles (alpha_f=0.38482)",y=1.02); fig.text(.02,-.02,note,fontsize=7); fig.tight_layout(); save(fig,"figure6_response_profiles")
print("matching plots written")
