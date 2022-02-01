#!/usr/bin/python3

import multiprocessing
import os
import modelcraft as mc
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


def _prepare_case(pdb_id):
    directory = os.path.join("data", "ep", pdb_id)
    if os.path.exists(directory):
        return None
    structure = pdbe.structure(pdb_id)
    rblocks = pdbe.rblocks(pdb_id)
    fmean, freer = sfdata.fmean_rfree(rblocks[0])
    phases = sfdata.phases(rblocks)
    if phases is None:
        return "No experimental phases deposited"
    if not sfdata.compatible_cell(structure, [fmean, freer, phases]):
        return "Different cell or space group in the structure and data"
    mc.update_cell(structure, new_cell=fmean.cell)
    refmac = mc.RefmacXray(structure, fmean, freer, cycles=10).run()
    if refmac.rfree > 0.06 * refmac.resolution_high + 0.17:
        return "R-free for deposited structure deemed too high"
    if refmac.data_completeness < 0.9:
        return "Data completeness less than 90%"
    phasematch = mc.PhaseMatch(fmean, phases, refmac.abcd).run()
    if phasematch.f_map_correlation < 0.2:
        return "F-map correlation less than 0.2"
    testset.write_case(pdb_id, directory, refmac, phasematch, fmean, freer, phases)


def prepare():
    pdb_ids = _search_for_pdb_ids()
    pdb_ids = {"1o6a", "2o7t", "4yme"}  # For small-scale testing
    pool = multiprocessing.Pool()
    failures = pool.map(_prepare_case, pdb_ids)
    testset.write_failures_table("prep_failures_ep.txt", failures)
