import gemmi
import tinterweb
from modelcraft.cell import remove_scale
from modelcraft.structure import (
    read_structure,
    remove_non_library_atoms,
    remove_residues,
)


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
    structure = read_structure(path)
    remove_residues(structure, "UNL")
    remove_non_library_atoms(structure)
    remove_scale(structure=structure)
    return structure
