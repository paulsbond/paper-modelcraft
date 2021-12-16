#!/usr/bin/python3

import argparse
import glob
import multiprocessing
import gemmi
import pandas
import requests
import solrq


_LOCK = multiprocessing.Lock()


def _request_json(url, data=None):
    _LOCK.acquire()
    if data is None:
        response = requests.get(url)
    else:
        response = requests.post(url, data=data)
    _LOCK.release()
    return response.json()


def _find_potential_pdbs():
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    query = solrq.Q(
        experimental_method="X-ray diffraction",
        experiment_data_available="y",
        resolution=solrq.Range(2, 2.3),
        number_of_polymer_entities=2,
        number_of_protein_chains=solrq.Range(1, solrq.ANY),
        number_of_DNA_chains=solrq.Range(1, solrq.ANY),
        uniprot_accession=solrq.Range(solrq.ANY, solrq.ANY),
        modified_residue_flag="N",
        max_observed_residues=solrq.Range(20, solrq.ANY),
    )
    filter_list = (
        "pdb_id,"
        "uniprot_accession,"
        "polymer_length,"
        "number_of_copies,"
        "overall_quality,"
    )
    request_data = {"q": query, "fl": filter_list, "rows": 1000000, "wt": "json"}
    response_json = _request_json(url, data=request_data)
    response_data = response_json.get("response", {})
    docs = response_data["docs"]
    for doc in docs:
        doc["uniprot_accession"] = doc["uniprot_accession"][0]
    data = pandas.DataFrame(docs)
    data.sort_values(by="overall_quality", ascending=False, inplace=True)
    data.drop_duplicates(subset="uniprot_accession", keep="first", inplace=True)
    return data.to_dict(orient="records")


def _check_alphafold_models(pdbs, db):
    for i in reversed(range(len(pdbs))):
        uniprot = pdbs[i]["uniprot_accession"]
        paths = glob.glob(f"{db}/AF-{uniprot}-*")
        if len(paths) != 1:
            del pdbs[i]
            continue
        doc = gemmi.cif.read(paths[0])
        block = doc.sole_block()


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("db", help="path to the AlphaFold human DB")
    args = parser.parse_args()
    pdbs = _find_potential_pdbs()
    print("Potential PDB entries:", len(pdbs))
    _check_alphafold_models(pdbs, args.db)
    print("With single AlphaFold models in the DB:", len(pdbs))


if __name__ == "__main__":
    _main()
