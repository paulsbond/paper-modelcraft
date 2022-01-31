#!/usr/bin/python3

import glob
import json
import multiprocessing
import os
import gemmi
import modelcraft as mc
import solrq
import pdbe
import tarfile
import tinterweb


def _make_search_structures(block, aligned):
    structure = gemmi.make_structure_from_block(block)
    assert len(structure) == 1
    model = structure[0]
    assert len(model) == 1
    chain = model[0]
    assert len(chain) == len(aligned)
    for i in reversed(range(len(chain))):
        if not aligned[i]:
            del chain[i]
    structures = {0: structure}
    table = block.find("_ma_qa_metric_local.", ["label_seq_id", "metric_value"])
    plddts = {int(row[0]): float(row[1]) for row in table}
    for cutoff in 50, 70, 90:
        structure = structure.clone()
        chain = structure[0][0]
        removed_any = False
        for i in reversed(range(len(chain))):
            plddt = plddts[chain[i].label_seq]
            if plddt < cutoff:
                del chain[i]
                removed_any = True
        if len(chain) == 0:
            break
        if removed_any:
            structures[cutoff] = structure
    return structures


def _get_sf_data(pdb_id):
    filename = f"r{pdb_id}sf.ent"
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    path = f"downloads/{filename}"
    _download_file(url, path)
    doc = gemmi.cif.read(path)
    rblocks = gemmi.as_refln_blocks(doc)
    cif2mtz = gemmi.CifToMtz()
    mtz = cif2mtz.convert_block_to_mtz(rblocks[0])
    ianom = next(DataItem.search(mtz, "KMKM"), None)
    imean = next(DataItem.search(mtz, "JQ"), None)
    fanom = next(DataItem.search(mtz, "GLGL"), None)
    fmean = next(DataItem.search(mtz, "FQ"), None)
    freer = FreeRFlag(ianom or imean or fanom or fmean).run().freer
    observations = ianom or imean or fanom or fmean
    return observations, freer


def _run_molrep(observations, structure, copies):
    molrep = Molrep(
        observations=observations,
        structure=structure,
        number_of_monomers=copies,
    ).run()
    return molrep


def _mean_plddt(structure):
    chain = structure[0][0]
    plddts = [atom.b_iso for residue in chain for atom in residue]
    return sum(plddts) / len(plddts)



# New


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


def _trim_alphafold(alphafold, residues):
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


def _superposed_completeness(trimmed, deposited, label_asym_id):
    ref_polymer = deposited[0].get_subchain(label_asym_id)
    wrk_polymer = trimmed[0].get_subchain("A")
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
    for pdb_id in pdb_ids & uniprot_data.keys():
        entry_data = uniprot_data[pdb_id]
        trimmed = _trim_alphafold(alphafold, entry_data["residues"])
        if trimmed is None:
            continue
        deposited = pdbe.structure(pdb_id)
        best_chain = entry_data["bestChainId"]
        completeness = _superposed_completeness(trimmed, deposited, best_chain)
        pdb_id = entry["accession"]
        if pdb_id in pdbs:
            mutated = False
            aligned = [False] * length
            for residues in entry["residues"]:
                if residues["indexType"] == "UNIPROT":
                    if "mutation" in residues:
                        mutated = True
                        break
                    start = residues["startIndex"] - 1
                    stop = residues["endIndex"]
                    for index in range(start, stop):
                        aligned[index] = True
            aligned = aligned[begin - 1 : end]
            if mutated or not any(aligned):
                del pdbs[pdb_id]
            else:
                pdbs[pdb_id]["aligned"] = aligned
    for key in list(pdbs.keys()):
        if "aligned" not in pdbs[key]:
            del pdbs[key]
    if len(pdbs) == 0:
        return None
    return max(pdbs.values(), key=lambda pdb: pdb["aligned"].count(True))


def _prepare_case(path):
    doc = gemmi.cif.read(path)
    block = doc.sole_block()
    uniprot = block.find_value("_ma_target_ref_db_details.db_accession")
    begin = int(block.find_value("_ma_target_ref_db_details.seq_db_align_begin"))
    end = int(block.find_value("_ma_target_ref_db_details.seq_db_align_end"))
    if begin != 1 or end == 1400:
        return
    try:
        alphafold = mc.read_structure(path)
    except (RuntimeError, ValueError):
        return
    pdb_ids, copies = _get_potential_pdbs(uniprot)
    if len(pdb_ids) == 0:
        return
    pdb_id = _choose_pdb(pdb_ids, uniprot, alphafold)
    # pdb_id = pdb["pdb_id"]
    # length = pdb["polymer_length"]
    # copies = pdb["number_of_copies"]
    # observations, freer = _get_sf_data(pdb_id)
    # metadata = {
    #     "uniprot": uniprot,
    #     "uniprot_begin": begin,
    #     "uniprot_end": end,
    #     "pdb_id": pdb_id,
    #     "pdb_length": length,
    #     "pdb_copies": copies,
    #     "pdb_resolution": observations.resolution_high(),
    # }
    # search_structures = _make_search_structures(block, pdb["aligned"])
    # best_mr_score = 0
    # best_structure = None
    # for cutoff, structure in search_structures.items():
    #     residues = len(structure[0][0])
    #     molrep = _run_molrep(observations, structure, copies)
    #     if molrep.n_solution == copies and molrep.mr_score > best_mr_score:
    #         best_mr_score = molrep.mr_score
    #         metadata["model_plddt_cutoff"] = cutoff
    #         metadata["model_plddt_mean"] = _mean_plddt(structure)
    #         metadata["model_residues"] = residues
    #         metadata["model_mr_score"] = molrep.mr_score
    #         metadata["model_mr_zscore"] = molrep.mr_zscore
    #         best_structure = molrep.structure
    # if best_mr_score == 0:
    #     return
    # contents = _entry_contents(pdb_id)
    # directory = f"data/{af_base}-{pdb_id}"
    # os.makedirs(directory, exist_ok=True)
    # contents.write_json_file(f"{directory}/contents.json")
    # write_mtz(f"{directory}/data.mtz", [observations, freer])
    # write_mmcif(f"{directory}/model.cif", best_structure)
    # with open(f"{directory}/metadata.json", "w") as stream:
    #     json.dump(metadata, stream, indent=4)

def prepare():
    paths = _alphafold_mmcif_paths()
    paths = [  # For small-scale testing
        "downloads/data/alphafold/AF-A4D1P6-F1-model_v2.cif.gz",
        "downloads/data/alphafold/AF-A6NK44-F1-model_v2.cif.gz",
        "downloads/data/alphafold/AF-P59666-F1-model_v2.cif.gz",
    ]
    pool = multiprocessing.Pool()
    pool.map(_prepare_case, paths)
