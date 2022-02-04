#!/usr/bin/python3

import glob
import json
import multiprocessing
import os
import shutil
import subprocess
import uuid
import gemmi
import modelcraft as mc
import pandas as pd
import pdbe
import ncacstat


_RESULTS_PATH = "results.csv"


def gather():
    done = None
    if os.path.exists(_RESULTS_PATH):
        results = pd.read_csv(_RESULTS_PATH)
        todo = results[results.isnull().any(axis=1)]
        done = results[~results.isnull().any(axis=1)]
        dirs = [f"data/{row['type']}/{row['id']}" for row in todo.to_records()]
    else:
        dirs = glob.glob("data/af/*") + glob.glob("data/ep/*") + glob.glob("data/mr/*")
    pool = multiprocessing.Pool()
    result_list = pool.map(_result, dirs)
    results = pd.DataFrame(result_list)
    if done is not None:
        results = pd.concat([done, results])
    results.to_csv(_RESULTS_PATH, index=False)


def _result(directory):
    split = directory.split("/")
    result = {"id": split[2], "type": split[1]}
    result.update(_metadata(directory))
    result.update(_ccp4i(directory))
    result.update(_modelcraft(directory))
    return result


def _metadata(directory):
    resolution = None
    f_map_correlation = None
    metadata_path = f"{directory}/metadata.json"
    if os.path.exists(metadata_path):
        with open(metadata_path) as stream:
            metadata = json.load(stream)
        resolution = metadata["data_resolution"]
        f_map_correlation = metadata["f_map_correlation"]
    return {"resolution": resolution, "f_map_correlation": f_map_correlation}


def _result_dict(key, completeness, rwork, rfree, seconds):
    return {
        f"{key}_completeness": completeness,
        f"{key}_rwork": rwork,
        f"{key}_rfree": rfree,
        f"{key}_seconds": seconds,
    }


def _ccp4i(directory):
    rwork = rfree = None
    seconds = 0
    log_path = f"{directory}/ccp4i/ccp4i.log"
    if os.path.exists(log_path):
        with open(log_path) as stream:
            for line in stream:
                if line[:12] == "Times: User:":
                    seconds += float(line.split()[2].rstrip("s"))
                if line[:19] == "           R factor":
                    rwork = float(line.split()[-1])
                if line[:19] == "             R free":
                    rfree = float(line.split()[-1])
    model_path = f"{directory}/ccp4i/ccp4i.pdb"
    completeness = _completeness(model_path)
    return _result_dict("ccp4i", completeness, rwork, rfree, seconds)


def _modelcraft(directory):
    rwork = rfree = seconds = None
    json_path = f"{directory}/modelcraft/modelcraft.json"
    if os.path.exists(json_path):
        with open(json_path) as stream:
            modelcraft = json.load(stream)
            if "final" in modelcraft and "termination_reason" in modelcraft:
                rwork = modelcraft["final"]["r_work"]
                rfree = modelcraft["final"]["r_free"]
                seconds = modelcraft["seconds"]["total"]
    model_path = f"{directory}/modelcraft/modelcraft.cif"
    completeness = _completeness(model_path)
    return _result_dict("modelcraft", completeness, rwork, rfree, seconds)


def _completeness(model_path):
    if not os.path.exists(model_path):
        return None
    pdb_id = model_path.split("/")[2]
    built = mc.read_structure(model_path)
    deposited = pdbe.structure(pdb_id)
    moved = _csymmatch(built, deposited)
    return ncacstat.completeness(moved, deposited)


def _csymmatch(structure, deposited):
    tmp_dir = f"tmp-{uuid.uuid4()}"
    os.mkdir(tmp_dir)
    ref_path = f"{tmp_dir}/ref.cif"
    wrk_path = f"{tmp_dir}/wrk.cif"
    out_path = f"{tmp_dir}/out.cif"
    _minimal_doc(deposited).write_file(ref_path)
    _minimal_doc(structure).write_file(wrk_path)
    args = ["csymmatch"]
    args += ["-pdbin-ref", ref_path]
    args += ["-pdbin", wrk_path]
    args += ["-pdbout", out_path]
    args += ["-origin-hand-work"]
    subprocess.call(args, stdout=subprocess.DEVNULL)
    moved = mc.read_structure(out_path)
    shutil.rmtree(tmp_dir)
    return moved


def _minimal_doc(structure):
    groups = gemmi.MmcifOutputGroups(False)
    groups.cell = True
    groups.symmetry = True
    groups.scale = True
    groups.atoms = True
    return structure.make_mmcif_document(groups)
