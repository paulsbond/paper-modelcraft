#!/usr/bin/python3

from typing import List, Optional
import distutils.spawn
import gemmi
import json
import math
import matplotlib.pyplot as plt
import multiprocessing
import numpy
import os
import pandas
import requests
import sklearn.linear_model
import solrq
import subprocess
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET


LOCK = multiprocessing.Lock()
INTERCEPT = 0.17
SLOPE = 0.06

class Entry:
    def __init__(self, pdbid: str):
        self.pdbid = pdbid
        self.resolution = None
        self.fmap = None
        self.rwork = None
        self.rfree = None
        self.failure_reason = None

    def path(self, filename: str) -> str:
        return os.path.join(self.pdbid, filename)

    def fail(self, reason: str) -> None:
        self.failure_reason = reason

    def download(self, filename: str) -> None:
        path = self.path(filename)
        if not os.path.exists(path):
            url = "https://www.ebi.ac.uk/pdbe/entry-files/download/" + filename
            LOCK.acquire()
            time.sleep(1.0)
            urllib.request.urlretrieve(url, path)
            LOCK.release()


def request_jcsg_entries() -> List[Entry]:
    request_url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    filter_list = "pdb_id,structure_determination_method"
    query = solrq.Q(
        experimental_method="X-ray diffraction",
        experiment_data_available="y",
        molecule_type="Protein",
        SG_center_name="JCSG",
    )
    request_data = {"fl": filter_list, "q": query, "rows": 1000000, "wt": "json"}
    print("Getting a list of SAD/MAD JCSG entries in the PDB")
    response = requests.post(request_url, data=request_data)
    if response.status_code != 200:
        print(f"Response with status code {response.status_code} received")
        sys.exit(response.text)
    response_data = response.json().get("response", {})
    pdbids = set()
    for doc in response_data["docs"]:
        if (
            "structure_determination_method" in doc
            and doc["structure_determination_method"] in (["SAD"], ["MAD"])
        ):
            pdbids.add(doc["pdb_id"])
    print(len(pdbids), "entries found")
    return [Entry(pdbid) for pdbid in pdbids]


def get_entries() -> List[Entry]:
    if os.path.exists("entries"):
        with open("entries") as stream:
            return [Entry(line.strip()) for line in stream]
    else:
        entries = request_jcsg_entries()
        with open("entries", "w") as stream:
            for entry in entries:
                stream.write(entry.pdbid + "\n")
        return entries


def make_freerflag_array(rblock: gemmi.ReflnBlock) -> numpy.ndarray:
    loop = rblock.block.find_loop("_refln.status")
    flags = []
    for status in loop:
        char = status[0]
        if char == "'" or char == '"':
            char = status[1]
        if char == "o":
            flags.append(1.0)
        elif char == "f":
            flags.append(0.0)
        else:
            flags.append(numpy.NAN)
    return numpy.array(flags, dtype=numpy.float32)


def block_is_compatible(mtz: gemmi.Mtz, rblock: gemmi.ReflnBlock) -> bool:
    if rblock.spacegroup is not None and rblock.spacegroup != mtz.spacegroup:
        return False
    if rblock.cell == gemmi.UnitCell():
        return True
    return (
        math.isclose(rblock.cell.a, mtz.cell.a, rel_tol=0.05)
        and math.isclose(rblock.cell.b, mtz.cell.b, rel_tol=0.05)
        and math.isclose(rblock.cell.c, mtz.cell.c, rel_tol=0.05)
        and math.isclose(rblock.cell.alpha, mtz.cell.alpha, rel_tol=0.05)
        and math.isclose(rblock.cell.beta, mtz.cell.beta, rel_tol=0.05)
        and math.isclose(rblock.cell.gamma, mtz.cell.gamma, rel_tol=0.05)
    )


def data_frame(entry: Entry, mtz: gemmi.Mtz, rblocks: gemmi.ReflnBlocks) -> Optional[pandas.DataFrame]:
    if "F_meas_au" not in rblocks[0].column_labels():
        return entry.fail("First block does not have F,SIGF")
    if not block_is_compatible(mtz, rblocks[0]):
        return entry.fail("First block has a different cell/spacegroup")
    hkl = rblocks[0].make_miller_array()
    df = pandas.DataFrame(data=hkl, columns=["H", "K", "L"])
    df["FREE"] = make_freerflag_array(rblocks[0])
    df["FP"] = rblocks[0].make_float_array("F_meas_au")
    df["SIGFP"] = rblocks[0].make_float_array("F_meas_sigma_au")
    abcd_rblocks = [b for b in rblocks if "pdbx_HL_A_iso" in b.column_labels()]
    if len(abcd_rblocks) == 0:
        return entry.fail("No ABCD blocks")
    abcd_rblock = abcd_rblocks[0]
    if not block_is_compatible(mtz, abcd_rblock):
        return entry.fail("ABCD block has a different cell/spacegroup")
    abcd_hkl = abcd_rblock.make_miller_array()
    abcd_df = pandas.DataFrame(data=abcd_hkl, columns=["H", "K", "L"])
    abcd_df["HLA"] = abcd_rblock.make_float_array("pdbx_HL_A_iso")
    abcd_df["HLB"] = abcd_rblock.make_float_array("pdbx_HL_B_iso")
    abcd_df["HLC"] = abcd_rblock.make_float_array("pdbx_HL_C_iso")
    abcd_df["HLD"] = abcd_rblock.make_float_array("pdbx_HL_D_iso")
    return df.merge(abcd_df, how="left", on=["H", "K", "L"])


def mtz_base(entry: Entry) -> gemmi.Mtz:
    cif_path = entry.path(f"{entry.pdbid}.cif")
    structure = gemmi.read_structure(cif_path)
    mtz = gemmi.Mtz()
    mtz.cell = structure.cell
    mtz.spacegroup = gemmi.SpaceGroup(structure.spacegroup_hm)
    mtz.add_dataset("HKL_base")
    return mtz


class RefmacResult:
    def __init__(self, path: str):
        xml = ET.parse(path).getroot()
        self.rwork = float(list(xml.iter("r_factor"))[-1].text)
        self.rfree = float(list(xml.iter("r_free"))[-1].text)
        self.resolution = float(list(xml.iter("resolution_high"))[-1].text)
        self.completeness = float(list(xml.iter("data_completeness"))[-1].text)


def rfree_is_too_high(entry: Entry) -> bool:
    cutoff = SLOPE * entry.resolution + INTERCEPT
    return entry.rfree > cutoff


def download_structure(entry: Entry) -> None:
    entry.download(f"{entry.pdbid}.cif")


def download_data(entry: Entry) -> None:
    entry.download(f"r{entry.pdbid}sf.ent")


def create_data_mtz(entry: Entry) -> None:
    mtz_path = entry.path("data.mtz")
    if os.path.exists(mtz_path):
        return
    mtz = mtz_base(entry)
    cif_path = entry.path(f"r{entry.pdbid}sf.ent")
    doc = gemmi.cif.read(cif_path)
    rblocks = gemmi.as_refln_blocks(doc)
    df = data_frame(entry, mtz, rblocks)
    if df is None:
        return
    column_labels = ["H", "K", "L", "FREE", "FP", "SIGFP", "HLA", "HLB", "HLC", "HLD"]
    column_types = ["H", "H", "H", "I", "F", "Q", "A", "A", "A", "A"]
    for column_label, column_type in zip(column_labels, column_types):
        mtz.add_column(column_label, column_type)
    mtz.set_data(df.to_numpy())
    mtz.update_reso()
    mtz.write_to_file(mtz_path)


def remove_unl(entry: Entry) -> None:
    xyzin = entry.path(f"{entry.pdbid}.cif")
    no_unl = entry.path("no_unl.cif")
    if os.path.exists(no_unl) or not os.path.exists(xyzin):
        return
    structure = gemmi.read_structure(xyzin)
    for model in structure:
        for chain in model:
            for i in reversed(range(len(chain))):
                if chain[i].name == "UNL":
                    del chain[i]
    structure.make_mmcif_document().write_file(no_unl)


def refmac(entry: Entry) -> None:
    hklin = entry.path("data.mtz")
    xyzin = entry.path("no_unl.cif")
    refmac_cif = entry.path("refmac.cif")
    refmac_err = entry.path("refmac.err")
    refmac_log = entry.path("refmac.log")
    refmac_mtz = entry.path("refmac.mtz")
    refmac_xml = entry.path("refmac.xml")
    if (
        os.path.exists(refmac_log)
        or not os.path.exists(hklin)
        or not os.path.exists(xyzin)
    ):
        return
    arguments = [
        "refmac5",
        "HKLIN", hklin,
        "XYZIN", xyzin,
        "HKLOUT", refmac_mtz,
        "XMLOUT", refmac_xml,
        "XYZOUT", refmac_cif,
    ]
    with open(refmac_log, "w") as log_stream:
        with open(refmac_err, "w") as err_stream:
            process = subprocess.Popen(
                args=arguments,
                stdin=subprocess.PIPE,
                stdout=log_stream,
                stderr=err_stream,
                encoding="utf8",
            )
            process.stdin.write("NCYCLES 10\n")
            process.stdin.write("MAKE NEWLIGAND NOEXIT\n")
            process.stdin.write("PHOUT\n")
            process.stdin.write("END\n")
            process.stdin.close()
            process.wait()


def check_refmac(entry: Entry) -> None:
    refmac_xml = entry.path("refmac.xml")
    if not os.path.exists(refmac_xml):
        return entry.fail("Refmac XML does not exist")
    result = RefmacResult(refmac_xml)
    if result.rfree == -999:
        return entry.fail("R-free value of -999")
    entry.resolution = result.resolution
    entry.rwork = result.rwork
    entry.rfree = result.rfree
    if rfree_is_too_high(entry):
        return entry.fail("High R-free value for the resolution")
    if result.completeness < 90:
        return entry.fail("Data completeness below 90%")
    refmac_cif = entry.path("refmac.cif")
    try:
        gemmi.read_structure(refmac_cif)
    except RuntimeError:
        return entry.fail("Refmac CIF file cannot be read by gemmi")


def cphasematch(entry: Entry) -> None:
    cphasematch_log = entry.path("cphasematch.log")
    cphasematch_err = entry.path("cphasematch.err")
    if os.path.exists(cphasematch_log):
        return
    tmp_mtz = entry.path("tmp.mtz")
    tmp_log = entry.path("tmp.log")
    tmp_err = entry.path("tmp.err")
    args = [
        "cmtzjoin",
        "-mtzout", tmp_mtz,
        "-mtzin", entry.path("data.mtz"),
        "-colin", "FP,SIGFP",
        "-colout", "FP,SIGFP",
        "-mtzin", entry.path("data.mtz"),
        "-colin", "HLA,HLB,HLC,HLD",
        "-colout", "HLA,HLB,HLC,HLD",
        "-mtzin", entry.path("refmac.mtz"),
        "-colin", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB",
        "-colout", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB",
    ]
    with open(tmp_log, "w") as log_stream:
        with open(tmp_err, "w") as err_stream:
            subprocess.call(args, stdout=log_stream, stderr=err_stream)
    args = [
        "cphasematch",
        "-mtzin", tmp_mtz,
        "-colin-fo", "FP,SIGFP",
        "-colin-hl-1", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB",
        "-colin-hl-2", "HLA,HLB,HLC,HLD",
    ]
    with open(cphasematch_log, "w") as log_stream:
        with open(cphasematch_err, "w") as err_stream:
            subprocess.call(args, stdout=log_stream, stderr=err_stream)
    if os.path.exists(tmp_mtz):
        os.remove(tmp_mtz)
    if os.path.exists(tmp_log):
        os.remove(tmp_log)
    if os.path.exists(tmp_err):
        os.remove(tmp_err)


def read_fmap_correlation(entry: Entry) -> None:
    cphasematch_log = entry.path("cphasematch.log")
    with open(cphasematch_log) as stream:
        lines = stream.readlines()
        for i, line in enumerate(lines):
            if line[:21] == "  Nrefl <fom1> <fom2>":
                split = lines[i+1].split()
                entry.fmap = float(split[6])
    if entry.fmap < 0.2:
        entry.fail("F-map correlation less than 0.2")


def get_sequence(entry: Entry) -> None:
    path = entry.path("sequence.fasta")
    if os.path.exists(path):
        return
    request_url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    filter_list = "entity_id,molecule_sequence"
    query = solrq.Q(pdb_id=entry.pdbid, molecule_type="Protein")
    request_data = {"fl": filter_list, "q": query, "rows": 1000000, "wt": "json"}
    LOCK.acquire()
    time.sleep(1)
    response = requests.post(request_url, data=request_data)
    LOCK.release()
    if response.status_code != 200:
        return entry.fail("No sequence data returned from the PDBe")
    response_data = response.json().get("response", {})
    with open(path, "w") as stream:
        for doc in response_data["docs"]:
            if (
                "entity_id" in doc
                and "molecule_sequence" in doc
                and len(doc["molecule_sequence"]) > 5
            ):
                stream.write(f">{doc['entity_id']}\n")
                for i in range(0, len(doc["molecule_sequence"]), 60):
                    stream.write(doc["molecule_sequence"][i : i + 60] + "\n")


def write_metadata(entries: List[Entry]) -> None:
    metadata = {}
    entries.sort(key=lambda entry: entry.pdbid)
    with open("failures", "w") as stream:
        for entry in entries:
            if entry.failure_reason is None:
                metadata[entry.pdbid] = {
                    "fmap": entry.fmap,
                    "resolution": entry.resolution,
                    "rfree": entry.rfree,
                    "rwork": entry.rwork,
                }
            else:
                line = f"{entry.pdbid} {entry.failure_reason}\n"
                stream.write(line)
    print("Fully processed", len(metadata), "entries")
    with open("metadata.json", "w") as stream:
        json.dump(metadata, stream, indent=2)


def prepare_entry(entry: Entry) -> Entry:
    os.makedirs(entry.pdbid, exist_ok=True)
    for step in (
        download_structure,
        download_data,
        create_data_mtz,
        remove_unl,
        refmac,
        check_refmac,
        cphasematch,
        read_fmap_correlation,
        get_sequence,
    ):
        step(entry)
        if entry.failure_reason is not None:
            return entry
    return entry


def prepare_entries(entries: List[Entry]) -> List[Entry]:
    for program in ("cmtzjoin", "cphasematch", "refmac5"):
        if not distutils.spawn.find_executable(program):
            sys.exit(f"Cannot find {program} executable")
    processes = os.cpu_count() or 1
    print("Preparing", len(entries), "entries using", processes, "processes")
    pool = multiprocessing.Pool(processes)
    return pool.map(prepare_entry, entries)


if __name__ == "__main__":
    entries = get_entries()
    entries = prepare_entries(entries)
    write_metadata(entries)
