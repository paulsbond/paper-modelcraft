#!/usr/bin/python3

import multiprocessing
import os
import modelcraft as mc
import environ
import pdbe
import sfdata
import testset


def _search_for_pdb_ids():
    query = (
        "molecule_type:Protein "
        "AND experimental_method:X\\-ray\\ diffraction "
        "AND experiment_data_available:y "
        "AND SG_center_name:JCSG "
        "AND ( "
        " structure_determination_method:SAD "
        " OR structure_determination_method:MAD "
        ") "
    )
    filter_list = "pdb_id"
    docs = pdbe.search(query, filter_list)
    return {doc["pdb_id"] for doc in docs}


def _fail(pdb_id, reason):
    testset.write_failure("ep", pdb_id, reason)


def _prepare_case(pdb_id):
    directory = os.path.join("data", "ep", pdb_id)
    if os.path.exists(directory) or testset.already_failed("ep", pdb_id):
        return
    print("Preparing", pdb_id)
    structure = pdbe.structure(pdb_id)
    rblocks = pdbe.rblocks(pdb_id)
    fmean, freer = sfdata.fmean_rfree(rblocks[0])
    phases = sfdata.phases(rblocks)
    if phases is None:
        return _fail(pdb_id, "No experimental phases deposited")
    if not sfdata.compatible_cell(structure, [fmean, freer, phases]):
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
    phasematch = mc.PhaseMatch(fmean, phases, refmac.abcd).run()
    if phasematch.f_map_correlation < 0.2:
        return _fail(pdb_id, "F-map correlation less than 0.2")
    testset.write_case(pdb_id, directory, refmac, phasematch, fmean, freer, phases)


def _prepare():
    environ.assert_ccp4()
    print("Preparing the EP testset...")
    pdb_ids = _search_for_pdb_ids()
    print(f"Found {len(pdb_ids)} potential entries")
    pool = multiprocessing.Pool()
    pool.map(_prepare_case, pdb_ids)
    testset.write_failures_table("ep")


if __name__ == "__main__":
    _prepare()
