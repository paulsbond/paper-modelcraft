#!/usr/bin/python3

import glob
import json
import os
import urllib.request
import gemmi
import pandas
import requests


_REQUEST_JSON_CACHE = {}


def _request_json(url, data=None):
    if not _REQUEST_JSON_CACHE and os.path.exists("request_cache.json"):
        with open("request_cache.json") as stream:
            _REQUEST_JSON_CACHE = json.load(stream)
    data = data or {}
    hash_ = hash((url, tuple(sorted(data.items()))))
    if hash_ in _REQUEST_JSON_CACHE:
        return _REQUEST_JSON_CACHE[hash_]
    if data is None:
        response = requests.get(url)
    else:
        response = requests.post(url, data=data)
    _REQUEST_JSON_CACHE[hash_] = response.json()
    with open("request_cache.json", "w") as stream:
        json.dump(_REQUEST_JSON_CACHE, stream)
    return response.json()


def _download(filename, url):
    path = os.path.join("downloads", filename)
    if not os.path.exists(path):
        os.makedirs("downloads", exist_ok=True)
        urllib.request.urlretrieve(url, path)
    return path


def _find_protein_na_pdbs(type_):
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    other_type = "RNA" if type_ == "DNA" else "DNA"
    query = (
        "experimental_method:X\\-ray\\ diffraction"
        " AND experiment_data_available:y"
        " AND resolution:[0 TO 4]"
        " AND number_of_polymer_entities:[2 TO 3]"  # allow for double stranded NA
        " AND number_of_protein_chains:[1 TO *]"
        f" AND number_of_{type_}_chains:[1 TO *]"
        f" AND -number_of_{other_type}_chains:[* TO *]"  # assure no other type
        " AND molecule_type:Protein"
        " AND uniprot_accession:[* TO *]"
        " AND modified_residue_flag:N"
        " AND max_observed_residues:[50 TO *]"
    )
    filter_list = "pdb_id,uniprot_accession"
    request_data = {"q": query, "fl": filter_list, "rows": 1000000, "wt": "json"}
    response_json = _request_json(url, data=request_data)
    response_data = response_json.get("response", {})
    docs = response_data["docs"]
    for doc in docs:
        doc["uniprot_accession"] = doc["uniprot_accession"][0]
    frame = pandas.DataFrame(docs)
    # Drop entries with two protein entities that match the above criteria
    frame.drop_duplicates(subset="pdb_id", keep=False, inplace=True)
    return list(frame.to_records(index=False))


def _uniprot_pdbs_dict(dna_pdbs, rna_pdbs):
    uniprot_entries = {}
    for pdb in dna_pdbs + rna_pdbs:
        pdb_id = pdb["pdb_id"]
        uniprot = pdb["uniprot_accession"]
        uniprot_entries.setdefault(uniprot, set()).add(pdb_id)
    return uniprot_entries


def _download_alphafold(uniprot):
    downloaded_models = glob.glob(os.path.join("downloads", f"AF-{uniprot}-*model*"))
    downloaded_errors = glob.glob(os.path.join("downloads", f"AF-{uniprot}-*error*"))
    if downloaded_models and downloaded_errors:
        return sorted(downloaded_models)[-1], sorted(downloaded_errors)[-1]
    server = "https://alphafold.ebi.ac.uk"
    url = f"{server}/api/search?type=main&start=0&rows=2"
    response_json = _request_json(url + "&q=uniprotAccession:" + uniprot)
    docs = response_json.get("docs", [])
    if len(docs) == 1:
        start = docs[0]["uniprotStart"]
        end = docs[0]["uniprotEnd"]
        if start == 1 and end != 1400:
            entry = docs[0]["entryId"]
            version = docs[0]["latestVersion"]
            model_filename = f"{entry}-model_v{version}.cif"
            error_filename = f"{entry}-predicted_aligned_error_v{version}.json"
            model_path = _download(model_filename, f"{server}/files/{model_filename}")
            error_path = _download(error_filename, f"{server}/files/{error_filename}")
            return model_path, error_path
    return None, None


def _pdbe_uniprot_data(uniprot):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/unipdb/{uniprot}"
    response_json = _request_json(url)
    data = response_json.get(uniprot, {}).get("data", [])
    return {entry["accession"]: entry for entry in data}


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


def _download_deposited(entry):
    filename = f"{entry}.cif"
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    return _download(filename, url)


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


def _create_test_case_for_uniprot(uniprot, pdbs):
    alphafold_model_path, alphafold_pae_path = _download_alphafold(uniprot)
    if alphafold_model_path is None or alphafold_pae_path is None:
        return
    try:
        alphafold = gemmi.read_structure(alphafold_model_path)
    except (RuntimeError, ValueError):
        return
    uniprot_data = _pdbe_uniprot_data(uniprot)
    completenesses = {}
    for entry in pdbs & uniprot_data.keys():
        entry_data = uniprot_data[entry]
        trimmed = _trim_alphafold(alphafold, entry_data["residues"])
        if trimmed is None:
            continue
        deposited_path = _download_deposited(entry)
        try:
            deposited = gemmi.read_structure(deposited_path)
        except (RuntimeError, ValueError):
            continue
        best_chain = entry_data["bestChainId"]
        completeness = _superposed_completeness(trimmed, deposited, best_chain)
        if 0.2 < completeness < 0.95:
            completenesses[entry] = completeness
    if completenesses:
        pdb = min(completenesses, key=completenesses.get)
        print("Chose PDB:", pdb)
        directory = os.path.join("data", "na", pdb)
        os.makedirs(directory, exist_ok=True)


def _main():
    dna_pdbs = _find_protein_na_pdbs("DNA")
    rna_pdbs = _find_protein_na_pdbs("RNA")
    print("Number of Protein/DNA PDB entries:", len(dna_pdbs))
    print("Number of Protein/RNA PDB entries:", len(rna_pdbs))
    uniprot_pdbs = _uniprot_pdbs_dict(dna_pdbs, rna_pdbs)
    print("Number of UniProt accessions:", len(uniprot_pdbs))
    for uniprot, pdbs in uniprot_pdbs.items():
        print("Processing UniProt:", uniprot)
        _create_test_case_for_uniprot(uniprot, pdbs)


if __name__ == "__main__":
    _main()
