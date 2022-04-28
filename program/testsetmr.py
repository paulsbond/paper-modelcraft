import multiprocessing
import os
import shutil
import tarfile
import urllib
import modelcraft as mc
import environ
import pdbe
import sfdata
import testset
import tinterweb


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


def _fail(key, reason):
    testset.write_failure("mr", key, reason)


def _prepare_case(pdb_id):
    directory = os.path.join("data", "mr", pdb_id)
    if os.path.exists(directory) or testset.already_failed("mr", pdb_id):
        return
    print("Preparing", pdb_id)
    try:
        structure = pdbe.structure(pdb_id)
    except urllib.error.HTTPError:
        return _fail(pdb_id, "Entry has been obsoleted in the PDB")
    rblocks = pdbe.rblocks(pdb_id)
    fmean, freer = sfdata.fmean_rfree(rblocks[0])
    if not sfdata.compatible_cell(structure, [fmean, freer]):
        return _fail(pdb_id, "Different cell or space group in the structure and data")
    mc.update_cell(structure, new_cell=fmean.cell)
    try:
        refmac = mc.RefmacXray(structure, fmean, freer, cycles=10).run()
    except ValueError:
        return _fail(pdb_id, "Refmac failure")
    if refmac.data_completeness < 90:
        return _fail(pdb_id, "Data completeness less than 90%")
    if refmac.rfree > 0.06 * refmac.resolution_high + 0.17:
        return _fail(pdb_id, "R-free for deposited structure deemed too high")
    model_path = f"downloads/data/reduced_full/{pdb_id.upper()}/model.pdb"
    model = mc.read_structure(model_path)
    model_refmac = mc.RefmacXray(model, fmean, freer, cycles=0).run()
    phasematch = mc.PhaseMatch(fmean, model_refmac.abcd, refmac.abcd).run()
    if phasematch.f_map_correlation < 0.2:
        return _fail(pdb_id, "F-map correlation less than 0.2")
    testset.write_case(pdb_id, directory, refmac, phasematch, fmean, freer)
    shutil.copy(model_path, f"{directory}/model.pdb")


def _prepare():
    environ.assert_ccp4()
    print("Preparing the MR testset...")
    pdb_ids = _pdb_ids()
    print(f"Found {len(pdb_ids)} potential entries")
    pool = multiprocessing.Pool()
    pool.map(_prepare_case, pdb_ids)
    testset.write_failures_table("mr")


if __name__ == "__main__":
    _prepare()
