import pandas as pd
from ase.io import read

data = pd.read_csv('solvent_prop.csv', names=['RC', 'a', 'b'])


for i in data.RC:
    try:
        A = read(i+'_extracted_rebuild.pdb')
        B = read(i+'_extracted_nosolv.pdb')
    except FileNotFoundError:
        continue
    # print(len(A), len(B))
    res = 'n'
    if len(A) > len(B):
        res = 'y'
    print(i, res)
