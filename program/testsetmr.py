import multiprocessing
import os
import shutil
import tarfile
import pdbe
import sfdata
import testset
import tinterweb
from modelcraft.cell import update_cell
from modelcraft.jobs.phasematch import PhaseMatch
from modelcraft.jobs.refmac import RefmacXray
from modelcraft.structure import read_structure


def _pdb_ids():
    directory = "downloads/data/reduced_full"
    if not os.path.exists(directory):
        server = "https://webfiles.york.ac.uk"
        filename = "Bond_et_al_2020_data.tar.gz"
        url = f"{server}/INFODATA/44145f0a-5d82-4604-9494-7cf71190bd82/{filename}"
        path = tinterweb.download_file(filename, url)
        tar = tarfile.open(path, "r:gz")
        tar.extract("data/full_reduced_set.tar.gz", path="downloads")
        tar.close()
        tar = tarfile.open("downloads/data/full_reduced_set.tar.gz", "r:gz")
        tar.extractall(path="downloads/data")
    return [pdb_id.lower() for pdb_id in os.listdir(directory)]


def _prepare_case(pdb_id):
    directory = os.path.join("data", "mr", pdb_id)
    if os.path.exists(directory):
        return None
    structure = pdbe.structure(pdb_id)
    rblocks = pdbe.rblocks(pdb_id)
    fmean, freer = sfdata.fmean_rfree(rblocks[0])
    if not sfdata.compatible_cell(structure, [fmean, freer]):
        return "Different cell or space group in the structure and data"
    update_cell(structure, new_cell=fmean.cell)
    refmac = RefmacXray(structure, fmean, freer, cycles=10).run()
    if refmac.rfree > 0.06 * refmac.resolution_high + 0.17:
        return "R-free deemed too high"
    if refmac.data_completeness < 0.9:
        return "Data completeness less than 90%"
    model_path = f"downloads/data/reduced_full/{pdb_id.upper()}/model.pdb"
    model = read_structure(model_path)
    model_refmac = RefmacXray(model, fmean, freer, cycles=0).run()
    phasematch = PhaseMatch(fmean, model_refmac.abcd, refmac.abcd).run()
    if phasematch.f_map_correlation < 0.2:
        return "F-map correlation less than 0.2"
    testset.write_case(pdb_id, directory, refmac, phasematch, fmean, freer)
    shutil.copy(model_path, f"{directory}/model.pdb")


def prepare():
    pdb_ids = _pdb_ids()
    pdb_ids = ["1bd9", "1bjn", "1e24"]  # For small-scale testing
    pool = multiprocessing.Pool()
    failures = pool.map(_prepare_case, pdb_ids)
    testset.write_failures_table("prep_failures_mr.txt", failures)
