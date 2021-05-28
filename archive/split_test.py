import stk
import glob


def split_molecule(mol, N, fg_end, core=False, fg='bromine'):
    """
    Split a molecule into N molecules and add functional group.

    Parameters
    ----------
    mol : :class:`stk.Molecule`
        Molecule to split.

    N : :class:`int`
        Number of molecules to split into. Each will contain at least
        one :attr:`fg_end` and :attr:`fg`.

    fg_end : :class:`str`
        Functional group to search for as starting point.

    fg : :class:`str`, optional
        Functional group to append at split point. Defaults to
        'bromine'.

    Returns
    -------
    molecules : :class:`list` of :class:`stk.Molecule`
        N molecules.

    """

    # TODO: Finish this function.
    molecules = []

    # Get number of fg_end.
    no_fg_end = 0
    if no_fg_end != N:
        raise ValueError(f'{N} {fg_end} were not found in molecule.')

    # For each fg_end, set a start atom.

    # Iterate through graph from

    if len(molecules) != N:
        raise ValueError(f'{N} molecules were not found.')

    return molecules


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
