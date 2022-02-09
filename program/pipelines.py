#!/usr/bin/python3

import glob
import multiprocessing
import os
import subprocess
import time
import gemmi
import modelcraft as mc


_LOCK = multiprocessing.Lock()


def _contents_path(directory):
    return os.path.join(directory, "contents.json")


def _sequence_path(directory):
    return os.path.join(directory, "sequence.fasta")


def _data_path(directory):
    return os.path.join(directory, "data.mtz")


def _model_path(directory):
    for filename in "model.pdb", "model.cif":
        path = os.path.join(directory, filename)
        if os.path.exists(path):
            return path
    return None


def _test_modelcraft(directory, disable=None):
    if disable is not None and directory.split("/")[1] != "mr":
        return
    modelcraft_dirname = f"modelcraft-no-{disable}" if disable else "modelcraft"
    modelcraft_dir = os.path.join(directory, modelcraft_dirname)
    if os.path.exists(modelcraft_dir):
        return
    args = ["modelcraft", "xray"]
    args += ["--directory", modelcraft_dir]
    args += ["--contents", _contents_path(directory)]
    args += ["--data", _data_path(directory)]
    if _model_path(directory):
        args += ["--model", _model_path(directory)]
    else:
        args += ["--unbiased"]
    if disable:
        args += [f"--disable-{disable}"]
    with _LOCK:
        time.sleep(1)
    subprocess.call(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _test_modelcraft_no_sheetbend(directory):
    _test_modelcraft(directory, disable="sheetbend")


def _test_modelcraft_no_pruning(directory):
    _test_modelcraft(directory, disable="pruning")


def _test_modelcraft_no_parrot(directory):
    _test_modelcraft(directory, disable="parrot")


def _test_modelcraft_no_dummy_atoms(directory):
    _test_modelcraft(directory, disable="dummy-atoms")


def _test_modelcraft_no_waters(directory):
    _test_modelcraft(directory, disable="waters")


def _test_modelcraft_no_side_chain_fixing(directory):
    _test_modelcraft(directory, disable="side-chain-fixing")


def _test_ccp4i(directory):
    ccp4i_dir = os.path.join(directory, "ccp4i")
    try:
        os.mkdir(ccp4i_dir)
    except FileExistsError:
        return
    mtz = gemmi.read_mtz_file(_data_path(directory))
    fsigf = next(mc.DataItem.search(mtz, "FQ"))
    freer = next(mc.DataItem.search(mtz, "I"))
    input_mtz = os.path.join(ccp4i_dir, "input.mtz")
    ccp4i_args = ["buccaneer_pipeline"]
    ccp4i_args += ["-seqin", _sequence_path(directory)]
    ccp4i_args += ["-mtzin", input_mtz]
    ccp4i_args += ["-colin-fo", fsigf.label()]
    ccp4i_args += ["-colin-free", freer.label()]
    ccp4i_args += ["-buccaneer-anisotropy-correction"]
    ccp4i_args += ["-buccaneer-fast"]
    ccp4i_args += ["-pdbout", os.path.join(ccp4i_dir, "ccp4i.pdb")]
    ccp4i_args += ["-prefix", os.path.join(ccp4i_dir, "ccp4i-")]
    if _model_path(directory):
        structure = mc.read_structure(_model_path(directory))
        refmac = mc.RefmacXray(structure, fsigf, freer, cycles=10).run()
        input_pdb = os.path.join(ccp4i_dir, "input.pdb")
        refmac.structure.write_pdb(input_pdb)
        mc.write_mtz(input_mtz, [fsigf, freer, refmac.abcd])
        ccp4i_args += ["-colin-hl", refmac.abcd.label()]
        ccp4i_args += ["-refmac-mlhl", "0"]
        ccp4i_args += ["-pdbin-mr", input_pdb]
        ccp4i_args += ["-buccaneer-keyword", "mr-model-seed"]
    else:
        contents = mc.AsuContents.from_json_file(_contents_path(directory))
        phases = next(mc.DataItem.search(mtz, "AAAA"))
        parrot = mc.Parrot(contents, fsigf, freer, phases).run()
        mc.write_mtz(input_mtz, [fsigf, freer, parrot.abcd])
        ccp4i_args += ["-colin-hl", parrot.abcd.label()]
        ccp4i_args += ["-refmac-mlhl", "1"]
    log_path = os.path.join(ccp4i_dir, "ccp4i.log")
    with open(log_path, "w") as log_stream:
        subprocess.call(ccp4i_args, stdout=log_stream, stderr=log_stream)


def run():
    print("Running pipelines...")
    dirs = glob.glob("data/af/*") + glob.glob("data/ep/*") + glob.glob("data/mr/*")
    dirs = ["data/mr/1bd9"]  # For small-scale testing
    pool = multiprocessing.Pool()
    pool.map_async(_test_modelcraft, dirs)
    pool.map_async(_test_modelcraft_no_sheetbend, dirs)
    pool.map_async(_test_modelcraft_no_pruning, dirs)
    pool.map_async(_test_modelcraft_no_parrot, dirs)
    pool.map_async(_test_modelcraft_no_dummy_atoms, dirs)
    pool.map_async(_test_modelcraft_no_waters, dirs)
    pool.map_async(_test_modelcraft_no_side_chain_fixing, dirs)
    pool.map_async(_test_ccp4i, dirs)
    pool.close()
    pool.join()
