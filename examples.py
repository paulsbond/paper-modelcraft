#!/usr/bin/python3

import argparse
import glob
import math
import os
import urllib.request
import gemmi
import pandas
import requests


def _request_json(url, data=None):
    if data is None:
        response = requests.get(url)
    else:
        response = requests.post(url, data=data)
    return response.json()


def _download(filename, url):
    path = os.path.join("downloads", filename)
    if not os.path.exists(path):
        os.makedirs("downloads", exist_ok=True)
        urllib.request.urlretrieve(url, path)
    return path


def _pdbe_uniprot_data(uniprot):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/unipdb/{uniprot}"
    response_json = _request_json(url)
    data = response_json.get(uniprot, {}).get("data", [])
    return {entry["accession"]: entry for entry in data}


def _download_deposited(entry):
    filename = f"{entry}.cif"
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    return _download(filename, url)


def _download_alphafold(uniprot):
    downloaded = glob.glob(os.path.join("downloads", f"AF-{uniprot}-*"))
    if downloaded:
        return downloaded[0]
    url = "https://alphafold.ebi.ac.uk/api/search?type=main&start=0&rows=2"
    response_json = _request_json(url + "&q=uniprotAccession:" + uniprot)
    docs = response_json.get("docs", [])
    if len(docs) == 1:
        start = docs[0]["uniprotStart"]
        end = docs[0]["uniprotEnd"]
        if start == 1 and end != 1400:
            entry = docs[0]["entryId"]
            version = docs[0]["latestVersion"]
            filename = f"{entry}-model_v{version}.cif"
            url = f"https://alphafold.ebi.ac.uk/files/{filename}"
            return _download(filename, url)
    return None


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


def _mean_plddt(structure):
    chain = structure[0][0]
    plddts = [atom.b_iso for residue in chain for atom in residue]
    return sum(plddts) / len(plddts)


def _find_protein_na_ids(type_, min_res, max_res):
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    other_type = "RNA" if type_ == "DNA" else "DNA"
    query = (
        "experimental_method:X\\-ray\\ diffraction"
        " AND experiment_data_available:y"
        f" AND resolution:[{min_res} TO {max_res}]"
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
    # Drop entries with two protein entities
    frame.drop_duplicates(subset="pdb_id", keep=False, inplace=True)
    return [tuple(row) for row in frame.to_numpy()]


def _superpose(trimmed, deposited, chain):
    trimmed.write_pdb("trimmed.pdb")
    deposited.write_pdb("deposited.pdb")
    # TODO: Use GESAMT
    return -1


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("min_res", type=float, help="minimum resolution")
    parser.add_argument("max_res", type=float, help="maximum resolution")
    args = parser.parse_args()
    dna_ids = _find_protein_na_ids("DNA", args.min_res, args.max_res)
    rna_ids = _find_protein_na_ids("RNA", args.min_res, args.max_res)
    uniprot_entries = {}
    for entry, uniprot in dna_ids + rna_ids:
        uniprot_entries.setdefault(uniprot, set()).add(entry)
    examples = {}
    for uniprot, entries in uniprot_entries.items():
        alphafold_path = _download_alphafold(uniprot)
        if alphafold_path is not None:
            alphafold = gemmi.read_structure(alphafold_path)
            uniprot_data = _pdbe_uniprot_data(uniprot)
            for entry in entries & uniprot_data.keys():
                entry_data = uniprot_data[entry]
                trimmed = _trim_alphafold(alphafold, entry_data["residues"])
                if trimmed is not None:
                    plddt = _mean_plddt(trimmed)
                    deposited_path = _download_deposited(entry)
                    deposited = gemmi.read_structure(deposited_path)
                    best_chain = entry_data["bestChainId"]
                    rmsd = _superpose(trimmed, deposited, best_chain)
                    example = {
                        "pdb": entry,
                        "chain": best_chain,
                        "uniprot": uniprot,
                        "plddt": round(plddt, 1),
                        "rmsd": round(rmsd, 3),
                    }
                    print(example)
                    break
            break


if __name__ == "__main__":
    _main()
