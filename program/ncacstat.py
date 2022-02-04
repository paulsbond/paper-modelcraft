import gemmi


_NAMES = {"N", "CA", "C"}
_TYPES = set()
_TYPES |= {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU"}
_TYPES |= {"LYS", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "UNK", "VAL"}


def completeness(
    built: gemmi.Structure, deposited: gemmi.Structure, cutoff: float = 1.0
) -> float:
    num_residues = 0
    num_residues_built = 0
    for ref_res in _protein_residues(deposited):
        num_residues += 1
        for wrk_res in _protein_residues(built):
            if all(_matching_atom(ref_res, wrk_res, name, cutoff) for name in _NAMES):
                num_residues_built += 1
                break
    return num_residues_built / num_residues


def _protein_residues(structure: gemmi.Structure):
    for chain in structure[0]:
        for residue in chain.first_conformer():
            if residue.name in _TYPES and all(name in residue for name in _NAMES):
                yield residue


def _matching_atom(res1, res2, name, cutoff):
    for atom1 in res1[name]:
        for atom2 in res2[name]:
            if atom1.pos.dist(atom2.pos) < cutoff:
                return True
    return False
