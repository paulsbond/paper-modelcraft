#!/usr/bin/python3

import glob
import json
import multiprocessing
import os
import shutil
import subprocess
import uuid
import gemmi
import pandas as pd


def _main():
    if "CCP4" not in os.environ:
        print("CCP4 environment variable not set")
        return
    dirs = glob.glob("data/af/*") + glob.glob("data/ep/*")
    pool = multiprocessing.Pool()
    results = pool.map(_result, dirs)
    pd.DataFrame(results).to_csv("results.csv", index=False)


def _result(directory):
    split = directory.split("/")
    result = {"id": split[2], "type": split[1]}
    result.update(_metadata(directory))
    result["ccp4i"] = _ccp4i_completeness(directory)
    result["modelcraft"] = _modelcraft_completeness(directory)
    return result


def _metadata(directory):
    metadata_path = f"{directory}/metadata.json"
    if os.path.exists(metadata_path):
        with open(metadata_path) as stream:
            return json.load(stream)
    return {"resolution": None, "f_map_correlation": None}


def _ccp4i_completeness(directory):
    return _completeness(directory, f"{directory}/ccp4i/ccp4i.pdb")


def _modelcraft_completeness(directory):
    return _completeness(directory, f"{directory}/modelcraft/modelcraft.cif")


def _completeness(directory, structure_path):
    deposited_path = f"{directory}/deposited.cif.gz"
    if os.path.exists(structure_path) or os.path.exists(deposited_path):
        try:
            structure = gemmi.read_structure(structure_path)
        except:
            return None
        try:
            deposited = gemmi.read_structure(deposited_path)
        except:
            return None
        return _csymmatch_ncacstat(structure, deposited)
    return None


def _csymmatch_ncacstat(structure, deposited):
    tmp_dir = f"tmp-{uuid.uuid4()}"
    os.mkdir(tmp_dir)
    ref_path = f"{tmp_dir}/ref.cif"
    deposited.make_mmcif_document().write_file(ref_path)
    wrk_path = f"{tmp_dir}/wrk.cif"
    structure.make_mmcif_document().write_file(wrk_path)
    args = ["csymmatch"]
    args += ["-pdbin-ref", "ref.cif"]
    args += ["-pdbin", "wrk.cif"]
    args += ["-pdbout", "sym.cif"]
    args += ["-origin-hand-work"]
    try:
        subprocess.call(args, cwd=tmp_dir, stdout=subprocess.DEVNULL)
    except:
        return None
    args = ["ncacstat", "ref.cif", "sym.cif"]
    try:
        output = subprocess.check_output(args, cwd=tmp_dir)
    except:
        return None
    shutil.rmtree(tmp_dir)
    split = output.split()
    if len(split) == 3:
        return float(split[1]) / float(split[0])
    return None


if __name__ == "__main__":
    _main()
