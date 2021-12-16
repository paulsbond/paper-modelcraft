#!/usr/bin/python3

import argparse
import glob
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


def _pdbe_uniprot_data(uniprot):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/unipdb/{uniprot}"
    response_json = _request_json(url)
    return response_json[uniprot]


def _download_mmcif(entry):
    filename = f"{entry}.cif"
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    path = os.path.join("downloads", filename)
    if not os.path.exists(path):
        os.makedirs("downloads", exist_ok=True)
        urllib.request.urlretrieve(url, path)
    return path


def _alphafold_structure(db_dir, uniprot):
    paths = glob.glob(os.path.join(db_dir, f"AF-{uniprot}-*"))
    if len(paths) != 1:
        return None
    doc = gemmi.cif.read(paths[0])
    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)
    return structure


def _mean_plddt(chain):
    plddts = [atom.b_iso for residue in chain for atom in residue]
    return sum(plddts) / len(plddts)


def _find_ids(type_, min_res, max_res):
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


def _align_alphafold_models(pdbs, db_dir):
    for i in reversed(range(len(pdbs))):
        uniprot = pdbs[i]["uniprot_accession"]
        paths = glob.glob(f"{db_dir}/AF-{uniprot}-*")
        if len(paths) != 1:
            del pdbs[i]
            continue
        doc = gemmi.cif.read(paths[0])
        block = doc.sole_block()
        begin = int(block.find_value("_ma_target_ref_db_details.seq_db_align_begin"))
        end = int(block.find_value("_ma_target_ref_db_details.seq_db_align_end"))
        aligned = _alignment(uniprot, pdbs[i]["pdb_id"], begin, end)
        if aligned is None or not any(aligned):
            del pdbs[i]
            continue
        structure = _make_search_structure(block, aligned)
        plddt = _mean_plddt(structure, block)
        pdbs[i]["aligned_plddt"] = plddt


def _alignment(uniprot, pdb, begin, end):
    length = response_json[uniprot]["length"]
    data = response_json[uniprot]["data"]
    entries = [e for e in data if e["accession"] == pdb]
    if len(entries) != 1:
        return None
    entry = entries[0]
    aligned = [False] * length
    for residues in entry["residues"]:
        if residues["indexType"] == "UNIPROT":
            if "mutation" in residues:
                return None
            start = residues["startIndex"] - 1
            stop = residues["endIndex"]
            for index in range(start, stop):
                aligned[index] = True
    return aligned[begin - 1 : end]


def _trim_chain(chain, aligned):
    assert len(chain) == len(aligned)
    for i in reversed(range(len(chain))):
        if not aligned[i]:
            del chain[i]
    return chain


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("db", help="path to a folder of AlphaFold models")
    parser.add_argument("min_res", type=float, help="minimum resolution")
    parser.add_argument("max_res", type=float, help="maximum resolution")
    args = parser.parse_args()
    dna_ids = _find_ids("DNA", args.min_res, args.max_res)
    rna_ids = _find_ids("RNA", args.min_res, args.max_res)
    dict_ = {}
    for entry, uniprot in dna_ids + rna_ids:
        dict_.setdefault(uniprot, []).append(entry)
    for uniprot, entries in dict_.items():
        af_structure = _alphafold_structure(args.db, uniprot)
        if af_chain is None:
            continue
        print(uniprot, entries)


if __name__ == "__main__":
    _main()
