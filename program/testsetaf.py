#!/usr/bin/python3

import glob
import multiprocessing
import os
import tarfile
import gemmi
import modelcraft as mc
import solrq
import environ
import pdbe
import sfdata
import testset
import tinterweb


def _alphafold_mmcif_paths():
    directory = "downloads/data/alphafold"
    if not os.path.exists(directory):
        server = "https://ftp.ebi.ac.uk"
        filename = "UP000005640_9606_HUMAN_v2.tar"
        url = f"{server}/pub/databases/alphafold/latest/{filename}"
        path = tinterweb.download_file(filename, url)
        tar = tarfile.open(path, "r")
        tar.extractall(path=directory)
    return glob.glob(f"{directory}/AF-*-F1-model_v2.cif.gz")


def _get_potential_pdbs(uniprot):
    query = solrq.Q(
        molecule_type="Protein",
        uniprot_accession=uniprot,
        modified_residue_flag="N",
        experimental_method="X-ray diffraction",
        experiment_data_available="y",
        resolution=solrq.Range(0, 4),
        number_of_polymer_entities=1,
        max_observed_residues=solrq.Range(20, solrq.ANY),
    )
    filter_list = "pdb_id,number_of_copies"
    docs = pdbe.search(query, filter_list)
    pdb_ids = {doc["pdb_id"] for doc in docs}
    copies = {doc["pdb_id"]: doc["number_of_copies"] for doc in docs}
    return pdb_ids, copies


def _truncate_alphafold(alphafold, residues):
    structure = alphafold.clone()
    chain = structure[0][0]
    observed = [False] * len(chain)
    for item in residues:
        if item.get("indexType") == "UNIPROT" and item.get("observed") == "Y":
            start = item["startIndex"] - 1
            stop = item["endIndex"]
            for i in range(start, stop):
                observed[i] = True
    if not any(observed):
        return None
    for i in reversed(range(len(chain))):
        if not observed[i]:
            del chain[i]
    return structure


def _superposed_similarity(truncated, deposited, label_asym_id):
    ref_polymer = deposited[0].get_subchain(label_asym_id)
    wrk_polymer = truncated[0].get_subchain("A")
    sup = gemmi.calculate_superposition(
        ref_polymer, wrk_polymer, gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP
    )
    ref_ca_positions = [
        residue["CA"][0].pos
        for residue in ref_polymer.first_conformer()
        if "CA" in residue
    ]
    sup_ca_positions = [
        gemmi.Position(sup.transform.apply(residue["CA"][0].pos))
        for residue in wrk_polymer.first_conformer()
        if "CA" in residue
    ]
    completed = 0
    for ref_pos in ref_ca_positions:
        for sup_pos in sup_ca_positions:
            if ref_pos.dist(sup_pos) < 1:
                completed += 1
                break
    return completed / len(ref_ca_positions)


def _choose_pdb(pdb_ids, uniprot, alphafold):
    uniprot_data = pdbe.uniprot_data(uniprot)
    chosen_pdb_id = None
    chosen_similarity = 0.9
    chosen_truncated = None
    chosen_deposited = None
    for pdb_id in pdb_ids & uniprot_data.keys():
        entry_data = uniprot_data[pdb_id]
        truncated = _truncate_alphafold(alphafold, entry_data["residues"])
        if truncated is None:
            continue
        deposited = pdbe.structure(pdb_id)
        best_chain = entry_data["bestChainId"]
        similarity = _superposed_similarity(truncated, deposited, best_chain)
        if 0.2 <= similarity <= chosen_similarity:
            chosen_pdb_id = pdb_id
            chosen_similarity = similarity
            chosen_truncated = truncated
            chosen_deposited = deposited
    return chosen_pdb_id, chosen_truncated, chosen_deposited, chosen_similarity


def _make_search_structures(structure):
    structures = [structure]
    for cutoff in 50, 70, 90:
        structure = structure.clone()
        chain = structure[0][0]
        removed_any = False
        for i in reversed(range(len(chain))):
            plddt = chain[i][0].b_iso
            if plddt < cutoff:
                del chain[i]
                removed_any = True
        if len(chain) == 0:
            break
        if removed_any:
            structures.append(structure)
    return structures


def _molecular_replacement(fmean, truncated, copies):
    best_mr_score = 0
    best_structure = None
    search_structures = _make_search_structures(truncated)
    for structure in search_structures:
        molrep = mc.Molrep(fmean, structure, copies).run()
        if molrep.n_solution == copies and molrep.mr_score > best_mr_score:
            best_mr_score = molrep.mr_score
            best_structure = molrep.structure
    return best_structure


def _fail(key, reason):
    testset.write_failure("af", key, reason)


def _prepare_case(path):
    doc = gemmi.cif.read(path)
    block = doc.sole_block()
    uniprot = block.find_value("_ma_target_ref_db_details.db_accession")
    directory = os.path.join("data", "af", uniprot)
    if os.path.exists(directory) or testset.already_failed("ep", uniprot):
        return
    print("Preparing", uniprot)
    end = int(block.find_value("_ma_target_ref_db_details.seq_db_align_end"))
    if end == 1400:
        return _fail(uniprot, "File contains residues 1-1400 of a longer sequence")
    alphafold = mc.read_structure(path)
    pdb_ids, copies = _get_potential_pdbs(uniprot)
    if len(pdb_ids) == 0:
        return _fail(uniprot, "No PDB entries")
    pdb_id, truncated, deposited, similarity = _choose_pdb(pdb_ids, uniprot, alphafold)
    if pdb_id is None:
        return _fail(uniprot, "No PDB entries with similarity between 20% and 90%")
    rblocks = pdbe.rblocks(pdb_id)
    fmean, freer = sfdata.fmean_rfree(rblocks[0])
    if not sfdata.compatible_cell(deposited, [fmean, freer]):
        return _fail(uniprot, "Different cell or space group in the structure and data")
    mc.update_cell(deposited, new_cell=fmean.cell)
    try:
        refmac = mc.RefmacXray(deposited, fmean, freer, cycles=10).run()
    except ValueError:
        return _fail(uniprot, "Refmac failure")
    if refmac.data_completeness < 0.9:
        return _fail(uniprot, "Data completeness less than 90%")
    if refmac.rfree > 0.06 * refmac.resolution_high + 0.17:
        return _fail(uniprot, "R-free for deposited structure deemed too high")
    mr_structure = _molecular_replacement(fmean, truncated, copies[pdb_id])
    if mr_structure is None:
        return _fail(uniprot, "Molecular replacement could not place all copies")
    molrep_refmac = mc.RefmacXray(mr_structure, fmean, freer, cycles=10).run()
    phasematch = mc.PhaseMatch(fmean, refmac.abcd, molrep_refmac.abcd).run()
    if phasematch.f_map_correlation < 0.2:
        return _fail(uniprot, "F-map correlation less than 0.2")
    testset.write_case(
        pdb_id,
        directory,
        refmac,
        phasematch,
        fmean,
        freer,
        superposed_similarity=similarity,
    )
    mc.write_mmcif(f"{directory}/model.cif", mr_structure)


def _prepare():
    environ.assert_ccp4()
    print("Preparing the AF testset...")
    paths = _alphafold_mmcif_paths()
    print(f"Found {len(paths)} potential entries")
    pool = multiprocessing.Pool()
    pool.map(_prepare_case, paths)
    testset.write_failures_table("af")


if __name__ == "__main__":
    _prepare()
