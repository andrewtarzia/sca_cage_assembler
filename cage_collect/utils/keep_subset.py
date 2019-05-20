import glob
import os
import random

N = 100


def check_prob(keep):
    problem_cifs = []
    if len(keep) == 0:
        return True
    for i in keep:
        if i in problem_cifs:
            return True
    return False


keep = []
CIFs = glob.glob('*.cif')
while check_prob(keep=keep):
    keep = random.sample(CIFs, N)

for i in CIFs:
    if i not in keep:
        # delete CIF
        os.remove(i)
