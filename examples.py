#!/usr/bin/python3

import argparse
import glob
import statistics
import gemmi
import pandas
import requests
import solrq


def _request_json(url, data=None):
    if data is None:
        response = requests.get(url)
    else:
        response = requests.post(url, data=data)
    return response.json()


def _find_potential_pdbs(type_, min_res, max_res):
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    query = solrq.Q(
        experimental_method="X-ray diffraction",
        experiment_data_available="y",
        resolution=solrq.Range(min_res, max_res),
        number_of_polymer_entities=2,
        number_of_protein_chains=solrq.Range(1, solrq.ANY),
        uniprot_accession=solrq.Range(solrq.ANY, solrq.ANY),
        modified_residue_flag="N",
        max_observed_residues=solrq.Range(20, solrq.ANY),
    )
    if type_ == "DNA":
        query &= solrq.Q(number_of_DNA_chains=solrq.Range(1, solrq.ANY))
    else:
        query &= solrq.Q(number_of_RNA_chains=solrq.Range(1, solrq.ANY))
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


def _align_alphafold_models(pdbs, db):
    for i in reversed(range(len(pdbs))):
        uniprot = pdbs[i]["uniprot_accession"]
        paths = glob.glob(f"{db}/AF-{uniprot}-*")
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
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/unipdb/{uniprot}"
    response_json = _request_json(url)
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


def _make_search_structure(block, aligned):
    structure = gemmi.make_structure_from_block(block)
    assert len(structure) == 1
    model = structure[0]
    assert len(model) == 1
    chain = model[0]
    assert len(chain) == len(aligned)
    for i in reversed(range(len(chain))):
        if not aligned[i]:
            del chain[i]
    return structure


def _mean_plddt(structure, block):
    table = block.find("_ma_qa_metric_local.", ["label_seq_id", "metric_value"])
    plddts = {int(row[0]): float(row[1]) for row in table}
    chain = structure[0][0]
    return statistics.mean(plddts[residue.label_seq] for residue in chain)


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("db", help="path to the AlphaFold human DB")
    parser.add_argument("type", choices=["DNA", "RNA"])
    parser.add_argument("min_res", type=float, help="minimum resolution")
    parser.add_argument("max_res", type=float, help="maximum resolution")
    args = parser.parse_args()
    pdbs = _find_potential_pdbs(args.type, args.min_res, args.max_res)
    print("Potential PDB entries:", len(pdbs))
    _align_alphafold_models(pdbs, args.db)
    print("Aligned with AlphaFold models in the DB:", len(pdbs))
    data = pandas.DataFrame(pdbs)
    data.sort_values(by="aligned_plddt", inplace=True)
    print(data)
    data.to_csv("examples.csv")


if __name__ == "__main__":
    _main()
