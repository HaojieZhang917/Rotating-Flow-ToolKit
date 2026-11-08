import numpy as np, sys
d=np.load(sys.argv[1]);
for k in ['x','u0','v0','w0','du0','dv0']: np.savetxt(sys.argv[2]+'_'+k+'.txt', d[k])
