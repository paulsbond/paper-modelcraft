import gemmi


_NAMES = {"N", "CA", "C"}


def completeness(built: gemmi.Structure, deposited: gemmi.Structure) -> float:
    num_residues = 0
    num_residues_built = 0
    search = gemmi.NeighborSearch(built, max_radius=1.0).populate(include_h=False)
    for chain in deposited[0]:
        for residue in chain.first_conformer():
            info = gemmi.find_tabulated_residue(residue.name)
            if info.kind == gemmi.ResidueInfoKind.AA:
                if all(name in residue for name in _NAMES):
                    num_residues += 1
                    if _residue_is_built(built, search, residue):
                        num_residues_built += 1
    return num_residues_built / num_residues


def _residue_is_built(built, search, residue):
    return all(_atom_is_built(built, search, residue, name) for name in _NAMES)


def _atom_is_built(built, search, residue, name):
    for atom_alt in residue[name]:
        for mark in search.find_atoms(atom_alt.pos, "\0", radius=1.0):
            cra = mark.to_cra(built[0])
            if cra.atom.name == name:
                return True
    return False
