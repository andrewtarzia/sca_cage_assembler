import pandas as pd
import logging
import sys
import os
from ase.visualize import view
from ase.io import read


def main():
    if (not len(sys.argv) == 4):
        print("""
    Usage: viz_structures_ase.py input_file done_file output_file
        input_file (str) - csv file to read in REFCODEs from
        done_file (str) - file with done REFCODEs
        output_file (str) - file to output selcted REFCODEs to
        """)
        sys.exit()
    else:
        input_file = sys.argv[1]
        done_file = sys.argv[2]
        output_file = sys.argv[3]

    done = []
    if os.path.isfile(done_file):
        with open(done_file, 'r') as f:
            for line in f.readlines():
                done.append(line.rstrip())
    else:
        open(done_file, 'x')

    to_keep = []
    if os.path.isfile(output_file):
        with open(output_file, 'r') as f:
            for line in f.readlines():
                to_keep.append(line.rstrip())
    else:
        open(output_file, 'x')

    logging.info(f'{len(done)} files done -- {len(to_keep)} to keep')

    data = pd.read_csv(input_file)
    RC = list(set(data.REFCODE))

    count = 0
    for i in sorted(RC):
        if i in done:
            count += 1
            continue
        logging.info('----------------------------')
        logging.info(f'{count} of {len(RC)}')
        RD = data[data.REFCODE == i]
        if RD.iloc[0].molecule == 0 and RD.iloc[0].pore_diam_opt == 0 and RD.iloc[0].no_windows == 0:
            done.append(i)
            count += 1
            continue
        logging.info(f'{RD}')
        max_pd = max(list(RD.pore_diam_opt))
        MD = RD[RD.pore_diam_opt == max_pd].iloc[0]
        logging.info(f'{MD.pore_diam_opt} -- {MD.no_windows}')
        file = str(MD.REFCODE)+'_MP_'+str(MD.molecule)+'_coms.pdb'
        logging.info(f'{file}')
        s = read(file)
        view(s)
        if input('keep? ') == 'w':
            to_keep.append(i)
    #    break
        count += 1
        done.append(i)
        logging.info('----------------------------')
        with open(done_file, 'w') as f:
            for l in done:
                f.write(l+'\n')
        with open(output_file, 'w') as f:
            for l in to_keep:
                f.write(l+'\n')


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='')
    main()
