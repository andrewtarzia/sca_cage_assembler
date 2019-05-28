import pandas as pd
from ase.visualize import view
from ase.io import read

done = []
with open('temp_done2.gcd', 'r') as f:
    for line in f.readlines():
        done.append(line.rstrip())

to_keep = []
with open('viz_select2.gcd', 'r') as f:
    for line in f.readlines():
        to_keep.append(line.rstrip())
print(len(done), len(to_keep))

data = pd.read_csv('classified.csv')
RC = list(set(data.REFCODE))

count = 0
#to_keep = []
for i in sorted(RC):
    if i in done:
        count += 1
        continue
    print('----------------------------')
    print(count, 'of', len(RC))
    RD = data[data.REFCODE == i]
    if RD.iloc[0].molecule == 0 and RD.iloc[0].pore_diam_opt == 0 and RD.iloc[0].no_windows == 0:
        done.append(i)
        count += 1
        continue
    print(RD)
    max_pd = max(list(RD.pore_diam_opt))
    MD = RD[RD.pore_diam_opt == max_pd].iloc[0]
    print(MD.pore_diam_opt, MD.no_windows)
    file = str(MD.REFCODE)+'_MP_'+str(MD.molecule)+'_coms.pdb'
    print(file)
    s = read(file)
    view(s)
    if input('keep? ') == 'w':
        to_keep.append(i)
#    break
    count += 1
    done.append(i)
    print('----------------------------')
    with open('temp_done2.gcd', 'w') as f:
        for l in done:
            f.write(l+'\n')
    with open('viz_select2.gcd', 'w') as f:
        for l in to_keep:
            f.write(l+'\n')

