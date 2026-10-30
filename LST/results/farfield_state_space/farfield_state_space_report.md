# Far-field first-order state-space validation

参数：R=285.36, alpha=0.38482, beta=0.07759, omega=-0.07759, N=69。

状态约定：X=[q;q_z]，q(z)=a exp(i*k_z*z)，lambda=i*k_z，E X_z=A X。

## Pencil count

- finite state modes: 6
- infinite generalized eigenvalues: 2
- original quadratic finite modes: 6
- finite modes bounded toward z=+infinity: 3
- left-null constraints of bounded finite subspace: 5
- direct left-null dimension of E: 1
- direct algebraic constraint rank: 1

## Root-by-root equivalence

|root|kz from state pencil|kz from quadratic pencil|difference|state residual|qz-lambda*q|quadratic residual|orientation under exp(i*kz*z)|
|---:|---:|---:|---:|---:|---:|---:|
|1|-0.0017175903295895709 - 0.39256794395253564im|-0.0017175903295895362 - 0.39256794395253547im|-3.469446951953614e-17 - 1.6653345369377348e-16im|2.937699515424621e-16|2.0659983184315618e-16|3.98407769182771e-17|grows toward +z; transfer toward wall under fixed edge normalization|
|2|8.480631267490856e-11 - 0.1491006262785122im|8.480638286449034e-11 - 0.1491006262785123im|-7.018958178471944e-17 + 1.1102230246251565e-16im|1.440602088915676e-16|1.4134633922413342e-16|4.491537174220476e-17|grows toward +z; transfer toward wall under fixed edge normalization|
|3|-2.3973648267336727e-10 - 0.14910062627753143im|-2.3973652543916577e-10 - 0.14910062627753104im|4.276579849744684e-17 - 3.885780586188048e-16im|2.652783311530037e-16|2.6001992701251985e-16|3.8537226635216804e-17|grows toward +z; transfer toward wall under fixed edge normalization|
|4|0.0017175903295824652 + 0.39256794395301686im|0.0017175903295825064 + 0.39256794395301686im|-4.119968255444917e-17 + 0.0im|4.628612174904732e-16|4.565544039024487e-16|7.883474335223003e-17|bounded toward +z|
|5|2.3974953789194126e-10 + 1.0335747363819936im|2.3974954222077454e-10 + 1.0335747363819863im|-4.328833272969723e-18 + 7.327471962526033e-15im|4.5152865374962875e-16|6.155005422945978e-16|2.490789001487803e-16|bounded toward +z|
|6|-8.48123979763955e-11 + 1.033574736383459im|-8.481239469627316e-11 + 1.0335747363834638im|-3.280122344490966e-18 - 4.6629367034256575e-15im|3.693307552620133e-16|5.140107649615208e-16|1.0191138891914541e-16|bounded toward +z|

## Direct sign check near the suspected root

The following values use the same `exp(i*k_z*z)` convention and evaluate the original 4×4 matrix directly.

|probe kz|minimum singular value of M_inf|det(M_inf)|
|---:|---:|---:|
|-0.0 - 0.149100626278im|1.0057133298625245e-13|4.605021159503932e-26 - 8.390112066459906e-28im|
|0.0 + 0.149100626278im|0.0006912381888134677|1.1265943699665275e-7 - 1.1520396558156657e-9im|

## Interpretation

The six finite roots from the first-order generalized pencil agree with the original quadratic roots to numerical precision. The two infinite eigenvalues are the algebraic constraints associated with the pressure-containing primitive formulation; they must not be silently discarded when a physical boundary operator is constructed.

Under q(z)=a exp(i*k_z*z), imag(k_z)>0 is bounded as z→+infinity and imag(k_z)<0 grows as z→+infinity. The latter roots can still describe a response propagated inward from a finite forcing plane, but they are not admissible outer-field decay modes by themselves.

This validation does not yet impose an incoming/outgoing scattering condition. The next implementation must select a physical finite/infinite constraint subspace and derive the correct number of independent outer boundary equations after the pressure constraint is treated.

## Bounded-subspace boundary-count warning

The current finite bounded subspace contains 3 modes, and its primitive-state left null space has 5 rows with residual 3.333088823870879e-16. Because the primitive descriptor pencil has 2 infinite eigenvalues but only 1 direct left-null algebraic constraint, at least one hidden higher-index constraint is present. These 5 rows cannot yet be inserted directly into the four-row collocation outer boundary block. A pressure/constraint reduction or an equivalent descriptor boundary formulation is required before solving the scattering problem.
