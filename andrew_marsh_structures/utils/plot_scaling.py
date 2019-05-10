import glob
from ase.io import read
import os
import matplotlib.pyplot as plt
succ = []
fail = []
for i in glob.glob('*_opt.xyz'):
    dir = i.replace('.xyz', '/')
    temp = dir + '.xtboptok'
    s = read(i)
    if os.path.isfile(temp):
        succ.append(len(s))
    else:
        fail.append(len(s))
# 3N*(3N+1)/2
xs = [(3/2)*i*(3*i+1) for i in succ]
xf = [(3/2)*i*(3*i+1) for i in fail]
plt.scatter(fail, xf, c='r')
plt.scatter(succ, xs, c='b')
plt.show()
