import pathlib
import json
import os


def read_json(file):
    with open(file, "r") as f:
        return json.load(f)


def read_envset_json(file):
    environment = read_json(file)
    for i in environment:
        if 'dir' in i or 'file' in i:
            environment[i] = pathlib.Path(environment[i])
            if 'dir' in i:
                if not os.path.exists(pathlib.Path(environment[i])):
                    os.mkdir(pathlib.Path(environment[i]))
    return environment


def xtb_path():
    return pathlib.Path(
        "/home/atarzia/miniconda3/envs/large_polyhedra/bin/xtb"
    )


def crest_path():
    raise NotImplementedError()
    return pathlib.Path("/home/atarzia/software/crest/crest")


def gulp_path():
    return pathlib.Path("/home/atarzia/software/gulp-6.1/Src/gulp")


def shape_path():
    return pathlib.Path(
        "/home/atarzia/software/shape_2.1_linux_64/"
        "SHAPE_2.1_linux_64/shape_2.1_linux64"
    )


def crest_conformer_settings(solvent=None):
    return {
        "conf_opt_level": "crude",
        "final_opt_level": "extreme",
        "charge": 0,
        "no_unpaired_e": 0,
        "max_runs": 1,
        "calc_hessian": False,
        "solvent": solvent,
        "crest_exec": crest_path(),
        "nc": 4,
        "etemp": 300,
        "keepdir": False,
        "cross": True,
        "md_len": None,
        "ewin": 5,
        "speed_setting": "squick",
    }
