#!/usr/bin/python3

import collections
import json
import multiprocessing
import os
import pdbe
import sfdata
from modelcraft.cell import update_cell
from modelcraft.jobs.refmac import RefmacXray
from modelcraft.jobs.phasematch import PhaseMatch
from modelcraft.scripts.contents import _entry_contents
from modelcraft.reflections import write_mtz


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
    os.makedirs(directory)
    structure = pdbe.structure(pdb_id)
    rblocks = pdbe.rblocks(pdb_id)
    fmean, freer = sfdata.fmean_rfree(rblocks[0])
    phases = sfdata.phases(rblocks)
    if phases is None:
        return "No experimental phases deposited"
    if not sfdata.compatible_cell(structure, [fmean, freer, phases]):
        return "Different cell or space group in the structure and data"
    update_cell(structure, new_cell=fmean.cell)
    refmac = RefmacXray(structure, fmean, freer, cycles=10).run()
    if refmac.rfree > 0.06 * refmac.resolution_high + 0.17:
        return "Rfree deemed too high"
    if refmac.data_completeness < 0.9:
        return "Data completeness less than 90%"
    phasematch = PhaseMatch(fmean, phases, refmac.abcd).run()
    if phasematch.f_map_correlation < 0.2:
        return "F-map correlation less than 0.2"
    with multiprocessing.Lock():
        contents = _entry_contents(pdb_id)
    metadata = {
        "data_resolution": refmac.resolution_high,
        "data_completeness": refmac.data_completeness,
        "deposited_rfree": refmac.rfree,
        "deposited_rwork": refmac.rwork,
        "f_map_correlation": phasematch.f_map_correlation,
    }
    data_mtz = os.path.join(directory, "data.mtz")
    contents_json = os.path.join(directory, "contents.json")
    sequence_fasta = os.path.join(directory, "sequence.fasta")
    metadata_json = os.path.join(directory, "metadata.json")
    write_mtz(data_mtz, fmean, freer, phases)
    contents.write_json_file(contents_json)
    contents.write_sequence_file(sequence_fasta)
    with open(metadata_json, "w") as stream:
        json.dump(metadata, stream, indent=2)


def prepare():
    pdb_ids = _search_for_pdb_ids()
    print("Found", len(pdb_ids), "possible entries for the PDB test set")
    pool = multiprocessing.Pool()
    failures = pool.map(_prepare_case, pdb_ids)
    counter = collections.Counter(failures)
    with open("ep.failures.txt", "w") as stream:
        for failure, count in counter.most_common():
            print(count, failure, file=stream)
