import gemmi
from modelcraft.cell import max_distortion
from modelcraft.jobs.ctruncate import CTruncate
from modelcraft.jobs.freerflag import FreeRFlag
from modelcraft.reflections import DataItem


def fmean_rfree(rblock):
    cif2mtz = gemmi.CifToMtz()
    mtz = cif2mtz.convert_block_to_mtz(rblock)
    ianom = next(DataItem.search(mtz, "KMKM"), None)
    imean = next(DataItem.search(mtz, "JQ"), None)
    fanom = next(DataItem.search(mtz, "GLGL"), None)
    fmean = next(DataItem.search(mtz, "FQ"), None)
    fmean = fmean or CTruncate(ianom or imean or fanom).run().fmean
    freer = DataItem(mtz, "FreeR_flag")
    freer_values = list(freer.columns[-1])
    freer_percentage = freer_values.count(0) / len(freer_values) * 100
    if freer_percentage == 0 or freer_percentage > 50:
        freer = FreeRFlag(fmean).run().freer
    return fmean, freer


def phases(rblocks):
    cif2mtz = gemmi.CifToMtz()
    for rblock in rblocks:
        mtz = cif2mtz.convert_block_to_mtz(rblock)
        abcd = next(DataItem.search(mtz, "AAAA"), None)
        if abcd is not None:
            return abcd


def compatible_cell(structure, mtzs):
    structure_spacegroup = gemmi.find_spacegroup_by_name(
        structure.spacegroup_hm,
        alpha=structure.cell.alpha,
        gamma=structure.cell.gamma,
    )
    for mtz in mtzs:
        if (
            structure_spacegroup.number != mtz.spacegroup.number
            or max_distortion(old_cell=structure.cell, new_cell=mtz.cell) > 0.05
        ):
            return False
    return True
