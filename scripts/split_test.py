import stk
import glob
from atools.stk_f import split_molecule


def main():

    Ns = [2, 3]
    for i, file in enumerate(glob.glob('*.mol')):
        print(file)
        bb = stk.BuildingBlock.init_from_file(
            file,
            functional_groups=['pyridine_N_metal']
        )
        print(bb.func_groups)
        splits = split_molecule(
            mol=bb,
            N=Ns[i],
            core=True,
            fg_end='pyridine_N_metal'
        )

        for j, spi in enumerate(splits):
            print(j)
            print(spi.func_groups)
            spi.write(f"{file.replace('.mol', '')}_{j}.mol")


if __name__ == '__main__':
    main()
