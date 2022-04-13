import gemmi
import modelcraft as mc


def fmean_rfree(rblock):
    cif2mtz = gemmi.CifToMtz()
    mtz = cif2mtz.convert_block_to_mtz(rblock)
    ianom = next(mc.DataItem.search(mtz, "KMKM"), None)
    imean = next(mc.DataItem.search(mtz, "JQ"), None)
    fanom = next(mc.DataItem.search(mtz, "GLGL"), None)
    fmean = next(mc.DataItem.search(mtz, "FQ"), None)
    if (ianom or imean or fanom or fmean) is None:
        raise RuntimeError("No observations found")
    fmean = fmean or mc.CTruncate(ianom or imean or fanom).run().fmean
    freer = mc.DataItem(mtz, "FreeR_flag")
    freer_values = list(freer.columns[-1])
    freer_percentage = freer_values.count(0) / len(freer_values) * 100
    if freer_percentage == 0 or freer_percentage > 50:
        freer = mc.FreeRFlag(fmean).run().freer
    return fmean, freer


def phases(rblocks):
    cif2mtz = gemmi.CifToMtz()
    for rblock in rblocks:
        try:
            mtz = cif2mtz.convert_block_to_mtz(rblock)
        except RuntimeError:  # Can occur if data not found in block
            continue
        abcd = next(mc.DataItem.search(mtz, "AAAA"), None)
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
            mtz.spacegroup is None  # e.g. 4zys
            or structure_spacegroup.number != mtz.spacegroup.number
            or mc.max_cell_distortion(old_cell=structure.cell, new_cell=mtz.cell) > 0.05
        ):
            return False
    return True
