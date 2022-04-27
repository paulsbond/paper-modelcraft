import gemmi
import modelcraft as mc
import tinterweb


def rblocks(pdb_id):
    filename = f"r{pdb_id}sf.ent"
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    path = tinterweb.download_file(filename, url)
    doc = gemmi.cif.read(path)
    return gemmi.as_refln_blocks(doc)


def search(query, filter_list):
    # Documentation: https://www.ebi.ac.uk/pdbe/api/doc/search.html
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    data = {"q": query, "fl": filter_list, "rows": 1000000, "wt": "json"}
    response = tinterweb.request_json(url, data)
    return response.get("response", {}).get("docs", [])


def structure(pdb_id):
    filename = f"{pdb_id}.cif"
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{filename}"
    path = tinterweb.download_file(filename, url)
    structure_ = mc.read_structure(path)
    mc.remove_residues(structure_, {"UNL", "UNX"})
    mc.remove_non_library_atoms(structure_)
    mc.remove_scale(structure=structure_)
    return structure_


def uniprot_data(uniprot):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/unipdb/{uniprot}"
    response_json = tinterweb.request_json(url)
    data = response_json.get(uniprot, {}).get("data", [])
    return {entry["accession"]: entry for entry in data}
