#!/usr/bin/python3

import glob
import multiprocessing
import os
import subprocess
import gemmi
from modelcraft.contents import AsuContents, PolymerType
from modelcraft.jobs.ctruncate import CTruncate
from modelcraft.jobs.parrot import Parrot
from modelcraft.jobs.refmac import RefmacXray
from modelcraft.reflections import DataItem, write_mtz
from modelcraft.structure import read_structure


def _contents_path(directory):
    return os.path.join(directory, "contents.json")


def _data_path(directory):
    return os.path.join(directory, "data.mtz")


def _model_path(directory):
    return os.path.join(directory, "model.cif")


def _fsigf(mtz):
    ianom = next(DataItem.search(mtz, "KMKM"), None)
    imean = next(DataItem.search(mtz, "JQ"), None)
    fanom = next(DataItem.search(mtz, "GLGL"), None)
    fmean = next(DataItem.search(mtz, "FQ"), None)
    return fmean or CTruncate(ianom or imean or fanom).run().fmean


def _test_modelcraft(directory):
    modelcraft_dir = os.path.join(directory, "modelcraft")
    try:
        os.mkdir(modelcraft_dir)
    except FileExistsError:
        return
    args = ["modelcraft", "xray"]
    args += ["--directory", modelcraft_dir]
    args += ["--contents", _contents_path(directory)]
    args += ["--data", _data_path(directory)]
    if os.path.exists(_model_path(directory)):
        args += ["--model", _model_path(directory)]
    else:
        args += ["--unbiased"]
    subprocess.call(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _test_ccp4i(directory):
    ccp4i_dir = os.path.join(directory, "ccp4i")
    try:
        os.mkdir(ccp4i_dir)
    except FileExistsError:
        return
    mtz = gemmi.read_mtz_file(_data_path(directory))
    fsigf = _fsigf(mtz)
    freer = next(DataItem.search(mtz, "I"), None)
    contents = AsuContents.from_file(_contents_path(directory))
    sequence_path = os.path.join(ccp4i_dir, "input.fasta")
    contents.write_sequence_file(sequence_path, types=[PolymerType.PROTEIN])
    input_mtz = os.path.join(ccp4i_dir, "input.mtz")
    ccp4i_args = ["buccaneer_pipeline"]
    ccp4i_args += ["-seqin", sequence_path]
    ccp4i_args += ["-mtzin", input_mtz]
    ccp4i_args += ["-colin-fo", fsigf.label()]
    ccp4i_args += ["-colin-free", freer.label()]
    ccp4i_args += ["-buccaneer-anisotropy-correction"]
    ccp4i_args += ["-buccaneer-fast"]
    ccp4i_args += ["-pdbout", os.path.join(ccp4i_dir, "ccp4i.pdb")]
    ccp4i_args += ["-prefix", os.path.join(ccp4i_dir, "ccp4i-")]
    if os.path.exists(_model_path(directory)):
        structure = read_structure(_model_path(directory))
        refmac = RefmacXray(structure, fsigf, freer, cycles=10).run()
        input_pdb = os.path.join(ccp4i_dir, "input.pdb")
        refmac.structure.write_pdb(input_pdb)
        write_mtz(input_mtz, [fsigf, freer, refmac.abcd])
        ccp4i_args += ["-colin-hl", refmac.abcd.label()]
        ccp4i_args += ["-refmac-mlhl", "0"]
        ccp4i_args += ["-pdbin-mr", input_pdb]
        ccp4i_args += ["-buccaneer-keyword", "mr-model-seed"]
    else:
        phases = DataItem(mtz, "HLA,HLB,HLC,HLD")
        parrot = Parrot(contents, fsigf, freer, phases).run()
        write_mtz(input_mtz, [fsigf, freer, parrot.abcd])
        ccp4i_args += ["-colin-hl", parrot.abcd.label()]
        ccp4i_args += ["-refmac-mlhl", "1"]
    subprocess.call(ccp4i_args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _main():
    if "CCP4" not in os.environ:
        print("CCP4 environment variable not set")
        return
    dirs = glob.glob("data/af/*") + glob.glob("data/ep/*")
    pool = multiprocessing.Pool()
    pool.map_async(_test_modelcraft, dirs)
    pool.map_async(_test_ccp4i, dirs)
    pool.close()
    pool.join()


if __name__ == "__main__":
    _main()
