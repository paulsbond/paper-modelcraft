import collections
import json
import multiprocessing
import os
import modelcraft as mc


_LOCK = multiprocessing.Lock()


def write_case(
    pdb_id,
    directory,
    refmac,
    phasematch,
    fmean,
    freer,
    phases=None,
    superposed_similarity=None,
):
    with _LOCK:
        contents = mc.entry_contents(pdb_id)
    metadata = {
        "data_resolution": refmac.resolution_high,
        "data_completeness": refmac.data_completeness,
        "deposited_rfree": refmac.rfree,
        "deposited_rwork": refmac.rwork,
        "f_map_correlation": phasematch.f_map_correlation,
    }
    if superposed_similarity is not None:
        metadata["superposed_similarity"] = superposed_similarity
    data_mtz = os.path.join(directory, "data.mtz")
    contents_json = os.path.join(directory, "contents.json")
    sequence_fasta = os.path.join(directory, "sequence.fasta")
    metadata_json = os.path.join(directory, "metadata.json")
    os.makedirs(directory)
    mc.write_mtz(data_mtz, [fmean, freer, phases])
    contents.write_json_file(contents_json)
    contents.write_sequence_file(sequence_fasta, types=[mc.PolymerType.PROTEIN])
    with open(metadata_json, "w") as stream:
        json.dump(metadata, stream, indent=4)


def write_failures_table(filename, failures):
    os.makedirs("tables", exist_ok=True)
    path = f"tables/{filename}"
    counter = collections.Counter(failures)
    with open(path, "w") as stream:
        for failure, count in counter.most_common():
            print(count, failure, file=stream)
