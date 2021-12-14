#!/usr/bin/python3

import json
import multiprocessing
import os
import subprocess
import gemmi
from modelcraft.contents import AsuContents, PolymerType
from modelcraft.jobs.parrot import Parrot
from modelcraft.jobs.refmac import RefmacXray
from modelcraft.reflections import DataItem, write_mtz
from modelcraft.structure import read_structure


def _paths(case):
    data = os.path.join("datasets", case, "data.mtz")
    contents = os.path.join("datasets", case, "contents.json")
    model = os.path.join("datasets", case, "model.pdb")
    return data, contents, model


def _call(args, log_path, err_path):
    with open(log_path, "w") as out_stream:
        with open(err_path, "w") as err_stream:
            subprocess.call(args=args, stdout=out_stream, stderr=err_stream)


def _test_ccp4i(case):
    directory = os.path.join("test", case, "ccp4i")
    os.makedirs(directory, exist_ok=True)
    data_path, contents_path, model_path = _paths(case)
    mtz = gemmi.read_mtz_file(data_path)
    fsigf = DataItem(mtz, "FP,SIGFP")
    freer = DataItem(mtz, "FREE")
    contents = AsuContents.from_file(contents_path)
    sequence_path = os.path.join(directory, "sequence.fasta")
    contents.write_sequence_file(sequence_path, types=[PolymerType.PROTEIN])
    ccp4i_dir = os.path.join(directory, "ccp4i")
    ccp4i_pdb = os.path.join(directory, "ccp4i.pdb")
    ccp4i_log = os.path.join(directory, "ccp4i.log")
    ccp4i_err = os.path.join(directory, "ccp4i.err")
    ccp4i_args = ["buccaneer_pipeline"]
    ccp4i_args += ["-seqin", sequence_path]
    ccp4i_args += ["-colin-fo", "FP,SIGFP"]
    ccp4i_args += ["-colin-free", "FREE"]
    ccp4i_args += ["-buccaneer-anisotropy-correction"]
    ccp4i_args += ["-buccaneer-fast"]
    ccp4i_args += ["-pdbout", ccp4i_pdb]
    ccp4i_args += ["-prefix", ccp4i_dir]
    if os.path.exists(model_path):
        structure = read_structure(model_path)
        refmac = RefmacXray(structure, fsigf, freer, cycles=10).run()
        refmac_pdb = os.path.join(directory, "refmac.pdb")
        refmac_mtz = os.path.join(directory, "refmac.mtz")
        refmac.structure.write_pdb(refmac_pdb)
        write_mtz(refmac_mtz, [fsigf, freer, refmac.abcd])
        ccp4i_args += ["-mtzin", refmac_mtz]
        ccp4i_args += ["-colin-hl", refmac.abcd.label()]
        ccp4i_args += ["-refmac-mlhl", "0"]
        ccp4i_args += ["-buccaneer-keyword", "mr-model-seed"]
        ccp4i_args += ["-pdbin-mr", refmac_pdb]
    else:
        phases = DataItem(mtz, "HLA,HLB,HLC,HLD")
        parrot = Parrot(contents, fsigf, freer, phases).run()
        parrot_mtz = os.path.join(directory, "parrot.mtz")
        write_mtz(parrot_mtz, [fsigf, freer, parrot.abcd])
        ccp4i_args += ["-mtzin", parrot_mtz]
        ccp4i_args += ["-colin-hl", parrot.abcd.label()]
        ccp4i_args += ["-refmac-mlhl", "1"]
    _call(ccp4i_args, ccp4i_log, ccp4i_err)


def _test_modelcraft(case):
    directory = os.path.join("test", case, "modelcraft")
    os.makedirs(directory, exist_ok=True)
    modelcraft_log = os.path.join(directory, "modelcraft.log")
    modelcraft_err = os.path.join(directory, "modelcraft.err")
    if os.path.exists(modelcraft_log):
        return
    data, contents, model = _paths(case)
    args = ["modelcraft", "xray"]
    args += ["--contents", contents]
    args += ["--data", data]
    if os.path.exists(model):
        args += ["--model", model]
    else:
        args += ["--unbiased"]
    args += ["--directory", directory]
    _call(args, modelcraft_log, modelcraft_err)


def _test(case):
    _test_ccp4i(case)
    # _test_modelcraft(case)


def _main():
    datasets = os.listdir("datasets")
    pool = multiprocessing.Pool()
    pool.map(_test, datasets)


if __name__ == "__main__":
    _main()
