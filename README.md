# Rotating Flow ToolKit

A Julia-based research toolkit for rotating-flow boundary layers, basic-flow
solvers, linear stability, and receptivity calculations.

## Packaged projects

The actively packaged rotating-disk workflow is located at:

```text
compress/BEK/Vonkarmen_bone
```

It provides Lopez generalized-Boussinesq and fully compressible stability
models, neutral-curve continuation, Sutherland marching, regression tests, and
usage documentation.

```bash
cd compress/BEK/Vonkarmen_bone
julia --project=.
```

See [`compress/BEK/Vonkarmen_bone/README.md`](compress/BEK/Vonkarmen_bone/README.md)
for the public API, scripts, tests, and numerical conventions.

Contact: hj_zhang@tju.edu.cn
