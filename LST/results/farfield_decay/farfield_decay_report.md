# Far-field decay analysis

配置：frozen-R finite-edge Dirichlet；R=285.36, alpha_f=0.38482, beta_f=0.07759, omega_f=-0.07759, N=69, Nk=16, kappa_max=12.0, zi=10.0, rf=0.0, c=0.5, n=0。

## Exponential fits

The fit is log(S)=const-gamma*zmax for amplitudes and log(E)=const-2*gamma*zmax for energies. No gamma was fixed during fitting.

|quantity|gamma|R2|RMS log residual|relative difference from 0.14910|
|---|---:|---:|---:|---:|
|S_f|-9.873220e-03|8.802846e-01|3.109107e-02|-1.066219e+00|
|Q_f|1.491072e-01|1.000000e+00|4.814219e-05|4.810254e-05|
|E_v|1.490339e-01|9.999999e-01|8.046557e-04|-4.434783e-04|
|E_w|-9.519271e-03|8.832988e-01|5.909230e-02|-1.063845e+00|
|E_full|1.456721e-01|9.996148e-01|4.883611e-02|-2.299071e-02|

Wall-gradient component amplitude fits:
component 1: gamma=-9.644189e-03, R2=8.826948e-01, relative difference=-1.064683e+00
component 2: gamma=-1.645551e-02, R2=8.194629e-01, relative difference=-1.110366e+00
component 3: gamma=6.297079e-02, R2=8.760297e-01, relative difference=-5.776607e-01

## Fixed-edge data

|zmax|S_f|Q_f|E_v|E_w|E_full|chi|
|---:|---:|---:|---:|---:|---:|---:|
|15.0|1.116286e-07|6.898565e-04|8.944607e-17|8.482779e-21|8.940353e-17|-5.704502e-04|
|20.0|1.227415e-07|3.272801e-04|2.020173e-17|1.013432e-20|2.018439e-17|-1.359257e-03|
|25.0|1.328889e-07|1.552881e-04|4.548093e-18|1.182680e-20|4.545687e-18|-3.121392e-03|
|30.0|1.394711e-07|7.368285e-05|1.023920e-18|1.300123e-20|1.029820e-18|-6.848371e-03|
|35.0|1.425498e-07|3.496213e-05|2.305673e-19|1.356783e-20|2.406946e-19|-1.409255e-02|
|40.0|1.427799e-07|1.658943e-05|5.198849e-20|1.360432e-20|6.396333e-20|-2.484225e-02|

## Single-kappa comparison

The Gaussian zmax sweep above is not the same problem as the earlier single-kappa=.4 sweep. The following fit isolates that difference using the same N=69 and zmax list.

|response|gamma|R2|RMS log residual|
|---|---:|---:|---:|
|E_v(single kappa)|1.492095e-01|9.999995e-01|1.821093e-03|
|E_w(single kappa)|-9.755907e-03|8.838155e-01|6.040938e-02|
|E_full(single kappa)|1.491410e-01|9.999990e-01|2.487509e-03|

## Far-field state fit

The primitive rows use [u,v,w,p]. The extended rows use [q,Dq] and modal derivatives i*k_z*a. The derivative scaling is a diagnostic first-order-state approximation; neither fit is treated as a rigorous physical projection. Near-degenerate roots are grouped as subspaces.

|zmax/state|fit residual|cond(V)|low vortical subspace|potential subspace|high vortical subspace|
|---|---:|---:|---:|---:|---:|
|15.0|1.121707e-02|5.083084e+01|9.960101e-01|2.352320e-06|3.987568e-03|
|extended q,Dq|2.468643e-02|2.437428e+01|9.946682e-01|2.687719e-06|5.329161e-03|
|20.0|1.121713e-02|5.083084e+01|9.960124e-01|4.662008e-07|3.987105e-03|
|extended q,Dq|2.468808e-02|2.437428e+01|9.946758e-01|5.878362e-07|5.323639e-03|
|25.0|1.121772e-02|5.083084e+01|9.960037e-01|2.440883e-07|3.996100e-03|
|extended q,Dq|2.470028e-02|2.437428e+01|9.946790e-01|3.949079e-07|5.320569e-03|
|30.0|1.122350e-02|5.083084e+01|9.960323e-01|7.058730e-07|3.967044e-03|
|extended q,Dq|2.475288e-02|2.437428e+01|9.947367e-01|2.005488e-06|5.261284e-03|
|35.0|1.125396e-02|5.083084e+01|9.960070e-01|4.183409e-07|3.992545e-03|
|extended q,Dq|2.501896e-02|2.437428e+01|9.948144e-01|6.171728e-06|5.179399e-03|
|40.0|1.139299e-02|5.083084e+01|9.958138e-01|1.404802e-05|4.172170e-03|
|extended q,Dq|2.619242e-02|2.437428e+01|9.947306e-01|4.731537e-06|5.264715e-03|

Far-field roots at the reference parameters:
root 1: kz=-2.3973652543916577e-10 - 0.14910062627753104im, type=vortical_like, pressure_ratio=4.039209e-12, vorticity_ratio=1.014196e+00, continuity_residual=1.474674e-17
root 2: kz=8.480638286449034e-11 - 0.1491006262785123im, type=vortical_like, pressure_ratio=4.347532e-17, vorticity_ratio=1.069699e+00, continuity_residual=5.657178e-18
root 3: kz=-0.0017175903295895362 - 0.39256794395253547im, type=potential_like, pressure_ratio=2.191754e-03, vorticity_ratio=2.792741e-09, continuity_residual=1.669756e-16
root 4: kz=0.0017175903295825064 + 0.39256794395301686im, type=potential_like, pressure_ratio=2.191680e-03, vorticity_ratio=2.792816e-09, continuity_residual=1.666780e-16
root 5: kz=2.3974954222077454e-10 + 1.0335747363819863im, type=vortical_like, pressure_ratio=3.792582e-13, vorticity_ratio=2.361655e+00, continuity_residual=5.632511e-17
root 6: kz=-8.481239469627316e-11 + 1.0335747363834638im, type=vortical_like, pressure_ratio=1.559975e-17, vorticity_ratio=2.816386e+00, continuity_residual=3.346544e-16

## Fixed-internal-Q diagnostic

The outer edge amplitude is complex-rescaled so that the local three-component velocity vector at z_ref=10 matches the zmax=20 reference in a least-squares sense. This is a diagnostic rescaling, not a physical incoming-wave boundary condition.

|zmax|complex edge scale|Q_f after scale|z_ref match error|wall shear|single-kappa Efull amplitude|
|---:|---:|---:|---:|---:|---:|
|15.0|0.47445320498071225 - 0.0003659682202343878im|3.273047e-04|8.822463e-04|5.296257e-08|2.187186e-08|
|20.0|1.0 + 0.0im|3.272801e-04|0.000000e+00|1.227415e-07|2.180515e-08|
|30.0|4.441642801295123 + 2.920643142206292e-5im|3.272729e-04|3.081518e-04|6.194808e-07|2.179875e-08|
|40.0|19.727833238757132 + 0.00012990840932749616im|3.272734e-04|2.258417e-03|2.816739e-06|2.184237e-08|

## Higher Chebyshev-order first-order check

This smaller high-order test fits qf only at zmax=20,30,40; it checks whether the fitted first-order rate changes with N.

|N|gamma_Q|gamma_S|relative gamma_Q difference from gamma0|
|---:|---:|---:|---:|
|69|1.491025e-01|-7.561218e-03|1.709833e-05|
|89|1.491027e-01|-1.102447e-02|1.781619e-05|
|109|1.491027e-01|-1.326297e-02|1.819912e-05|

## Exact kappa=0 pressure gauge check

The original direct solve rejects the exact zero mode. An explicit pseudoinverse gauge was tested only after diagnosing the singular vector.

|quantity|value|
|---|---:|
|smallest singular value|1.371233e-14|
|second smallest singular value|1.931242e-04|
|left-null compatibility ratio|1.384747e-12|
|gauged interior residual|1.107386e-12|
|gauged boundary residual|1.118016e-14|
|gauged velocity norm|1.380799e+00|
|right-null velocity norm|2.285249e-13|
|right-null pressure norm|1.000000e+00|

## Interpretation

Answer 1: the internal first-order norm Q_f is controlled almost entirely by the two near-degenerate low-|imag(kz)| vortical roots kz≈-0.14910i in the primitive-state diagnostic. They are transfer roots under the fixed-edge normalization, not bounded outer-field decay modes. The fit residual decreases as the outer window moves outward, but the modal matrix becomes ill-conditioned; this supports subspace compatibility rather than a rigorous projection.
Answer 2: the single-kappa=.4 full coupling energy has amplitude-equivalent gamma_full≈0.14914, matching the low-|imag(kz)| transfer rate. The Gaussian-spectrum full energy has gamma_full≈0.14567 because the wall contribution and finite spectral reconstruction do not follow the same single exponential.
Answer 3: the volume contribution follows gamma_v≈0.14903, while the wall contribution has gamma_w<0 in the fit and approaches a plateau. The two mechanisms therefore do not share one universal decay law.
Answer 4: fixed-edge zmax changes both the distance from the Dirichlet forcing to the wall and the finite-domain boundary-value problem. After rescaling the edge input to keep the internal Q_f norm fixed, the single-kappa volume response changes by less than about 0.3%, whereas the wall response still grows strongly.
Answer 5: gamma_Q is stable from N=69 to N=109, but gamma_S changes with N and the wall-gradient magnitude is not converged. The internal/volume decay conclusion is robust at this stage; the wall and total Gaussian response require additional wall-gradient and Fourier convergence.
## Sign and normalization audit

The far-field convention is q(z)=a*exp(i*kz*z) with the same exp(i*alpha*r+i*beta*theta-i*omega*t) convention as the spatial operator. Therefore |exp(i*kz*z)|=exp(-imag(kz)*z): imag(kz)>0 decays toward z=+infinity, while imag(kz)<0 grows toward z=+infinity.

The fitted root kz≈-0.14910i must therefore not be called a bounded outer-field decay mode. In the fixed-edge calculation it is normalized at z=zmax and propagated toward the wall: q(z)=q(zmax)*exp(0.14910*(z-zmax)), so the wall/interior amplitude decreases as zmax is moved outward. This is a transfer decay with forcing distance, not yet an incoming-wave boundary condition at infinity.

For amplitudes the fit is log(A)=const-gamma*zmax. For quadratic energies the code fits log(E)=const-2*gamma*zmax and reports gamma=-slope/2. Thus gamma_Q, gamma_Ev and gamma_full are compared with the far-field root on the same amplitude-decay convention.

The primitive-state fit is not a rigorous first-order-system projection because the current far-field API returns only [u,v,w,p] eigenvectors, not the complete derivative state. Near-degenerate modes must be interpreted as subspaces.
The coupled response is reconstructed with physical epsilon and the fixed edge amplitude already contained in qf; no additional a or epsilon factor is applied after reconstruction.
