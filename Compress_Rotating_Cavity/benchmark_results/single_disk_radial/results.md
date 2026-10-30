# Single-disk radial attraction experiment

Steady axisymmetric boundary-layer model; all cases retain the radial non-similar derivatives.

Configuration: degree=80, eta_max=20.0, actual ds=0.01993580167094412, requested r/r0=100.00000000000004, Newton residual tolerance=2.0e-10.

| Inlet | Status | Initial relative FG error | Final relative FG error | Fitted spatial exponent |
|---|---|---:|---:|---:|
| similarity | completed | 5.797489e-24 | 2.430829e-16 | NaN |
| small_shape | completed | 9.560122e-03 | 6.551543e-11 | -3.83617 |
| small_opposite | completed | 4.478379e-03 | 6.782262e-11 | -3.84707 |
| large_radial | completed | 3.902578e-01 | 5.717870e-09 | -3.92560 |
| large_swirl | completed | 2.160691e-01 | 9.954875e-10 | -3.82386 |

Fitted exponents use 0.8 <= log(r/r0) <= 1.8; they are finite-interval fits, not eigenvalues. Values below max(100 times the control drift, 1e-10) are excluded.
Late-stage error plateaus can reflect the nonlinear stopping tolerance. Errors are relative to the reference on the same finite eta grid; they are not absolute errors relative to the infinite-domain solution.
The inlet has prescribed F,G only. Its H entries are NaN because inlet H is not independently prescribed or inferred by imposing F_s=0.
For s>0, H is solved with H_eta+2F+F_s=0. H at the far boundary is not prescribed.

These finite inlet experiments test attraction within the steady axisymmetric outward-flow boundary-layer model. They do not prove global attraction for arbitrary inlets, finite disk NS flows, or three-dimensional physical stability.
The r/r0 interval is a mathematical continuation interval; no physical local Reynolds-number interval has been assigned.
No reverse velocity clipping, FLARE, or projection onto a similarity profile is used. A converged reverse flow terminates the run explicitly.

Verification results are written separately by scripts/verify_single_disk_radial.jl.
