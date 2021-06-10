
def xtb_path():
    return '/home/atarzia/anaconda3/envs/sca_building/bin/xtb'


def crest_path():
    return '/home/atarzia/software/crest/crest'


def gulp_path():
    return '/home/atarzia/software/gulp-5.1/Src/gulp/gulp'


def shape_path():
    return (
        '/home/atarzia/software/shape_2.1_linux_64/'
        'SHAPE_2.1_linux_64/shape_2.1_linux64'
    )


def crest_conformer_settings(solvent=None):
    return {
        'conf_opt_level': 'crude',
        'final_opt_level': 'extreme',
        'charge': 0,
        'no_unpaired_e': 0,
        'max_runs': 1,
        'calc_hessian': False,
        'solvent': solvent,
        'crest_exec': crest_path(),
        'nc': 4,
        'etemp': 300,
        'keepdir': False,
        'cross': True,
        'md_len': None,
        'ewin': 5,
        'speed_setting': 'squick',
    }
