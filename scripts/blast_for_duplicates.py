#!/usr/bin/env python3

from argparse import ArgumentParser, FileType, SUPPRESS
import glob
import json
import os
import shutil
import subprocess
import tempfile
from typing import Dict, List, Tuple

from Bio import SeqIO

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXCEPTIONS = (
    {"BGC0000214", "BGC0001074"},  # S. olindensis DAUPE 5622 / DAUPE 5622 2
    {"BGC0000849", "BGC0000038"},  # 849 is contained in 38, but a different compound
    {"BGC0002040", "BGC0002462"},  # Different strains
    {"BGC0002040", "BGC0002282"},  # Different strains
    {"BGC0002407", "BGC0002621"},  # Different strains
    {"BGC0001011", "BGC0001013"},  # Different strains
    {"BGC0000100", "BGC0001670"},  # Different strains
    {"BGC0001023", "BGC0001429"},  # Different strains
    {"BGC0002294", "BGC0002313"},  # Different strains
    {"BGC0001956", "BGC0002479"},  # Different strains
    {"BGC0000703", "BGC0000705"},  # Different strains
    {"BGC0000720", "BGC0000721"},  # Different strains
    {"BGC0000742", "BGC0000743"},  # Different strains
    {"BGC0000744", "BGC0000745"},  # Different strains
    {"BGC0000746", "BGC0000748"},  # Different strains
    {"BGC0000693", "BGC0000694"},  # Possibly different strains
    {"BGC0000416", "BGC0001692"},  # Clusters overlap
)
ALREADY_RETIRED = (
    {"BGC0000246", "BGC0000247"},
    {"BGC0002225", "BGC0002728"},
    {"BGC0001539", "BGC0001673"},
    {"BGC0000672", "BGC0002388"},
    {"BGC0001994", "BGC0002094"},
    {"BGC0000943", "BGC0002463"},
    {"BGC0000821", "BGC0000822"},
    {"BGC0000821", "BGC0000823"},
    {"BGC0000822", "BGC0000823"},
    {"BGC0001985", "BGC0002041"},
    {"BGC0001985", "BGC0002292"},
    {"BGC0002041", "BGC0002292"},
    {"BGC0000270", "BGC0000835"},
)


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("files", nargs="*",
                        help="File(s) to check for duplicate sequences")
    parser.add_argument("-c", "--cutoff", type=float, default=97.0,
                        help="Similarity cutoff to report, in percent")
    parser.add_argument("-l", "--length-difference-cutoff", type=float, default=0.05,
                        help="Length difference cutoff to use (default: %(default)s)")
    parser.add_argument("-d", "--debug", action="store_true", default=False,
                        help="Save debug info")
    parser.add_argument("-a", "--analyse-only", type=FileType("r", encoding="utf-8"), default=SUPPRESS,
                        help="Skip running blast and just analyse an existing output table")
    args = parser.parse_args()

    if not args.files:
        args.files = sorted(glob.glob(os.path.join(BASE_DIR, "data", "*.json")))

    if "analyse_only" in args:
        analyse_blast(args.analyse_only.read(), args.cutoff, args.length_difference_cutoff)
    else:
        run(args.files, args.cutoff, args.length_difference_cutoff, args.debug)


IDMap = Dict[str, Tuple[str, int, int]]


def run(infiles: List[str], cutoff: float, length_diff: float, debug: bool = True) -> None:
    id_map = extract_ids(infiles)
    with tempfile.NamedTemporaryFile(suffix=".db") as db_handle:
        db_name = db_handle.name
        fasta_map = create_fastas(id_map)
        fastas = "\n".join(fasta_map.values())

        if debug:
            with open("all_seqs.fa", "w", encoding="utf-8") as handle:
                handle.write(fastas)

        cmd = ["makeblastdb", "-dbtype", "nucl", "-title", "mibig_validation", "-out", db_name]
        result = subprocess.run(cmd, input=fastas, capture_output=True, encoding="utf-8")
        result.check_returncode()

        if debug:
            for fname in glob.glob(f"{db_name}.*"):
                _, ext = os.path.splitext(fname)
                shutil.copyfile(fname, f"all_seqs{ext}")
        
        cmd = ["blastn", "-db", db_name, "-outfmt", "6 qacc sacc qlen slen length pident", "-num_threads", "8"]
        result = subprocess.run(cmd, input=fastas, capture_output=True, encoding="utf-8")
        result.check_returncode()

        if debug:
            with open("all_seqs.tbl", "w", encoding="utf-8") as handle:
                handle.write(result.stdout)

        analyse_blast(result.stdout, cutoff, length_diff)


def analyse_blast(result: str, cutoff: float, length_diff: float):
    for line in result.split("\n"):
        if not line:
            break
        qacc, sacc, qlen_raw, slen_raw, length_raw, pident_raw = line.split("\t")
        if qacc == sacc:
            continue
        qlen = int(qlen_raw)
        slen = int(slen_raw)
        length = int(length_raw)
        pident = float(pident_raw)
        diff = abs(length - qlen)
        if pident >= cutoff and diff / qlen <= length_diff:
            qmacc, qracc = qacc.split("|")[:2]
            smacc, sracc = sacc.split("|")[:2]
            qdata = get_entry_data(qmacc)
            sdata = get_entry_data(smacc)
            shared_pubs = sorted(qdata.publications & sdata.publications)
            shared_compounds = sorted(qdata.compounds & sdata.compounds)
            if {qmacc, smacc} in EXCEPTIONS or {qmacc, smacc} in ALREADY_RETIRED:
                continue
            if shared_pubs or shared_compounds or qdata.name == sdata.name:
                prefix = get_common_prefix([qdata.name, sdata.name])
                qunique = qdata.name[len(prefix):]
                sunique = sdata.name[len(prefix):]
                name = prefix
                if not qunique and not sunique:
                    print(f"{qmacc} <-> {smacc}: {pident}\t{qracc:>14}:{qlen:08d}\t{sracc:>14}:{slen:08d}\t{diff}\t{round(diff / qlen, 2)}\t{name}\t{shared_pubs}\t{shared_compounds}")


class EntryData:
    def __init__(self, publications: List[str], compounds: List[str], organism_name: str) -> None:
        self.publications = set(publications)
        self.compounds = set(compounds)
        self.name = organism_name

    def __str__(self) -> str:
        return ", ".join(self.publications) + " : " + ", ".join(self.compounds)


def get_common_prefix(names: List[str]):
    if not names:
        return ""
    shortest = min(names)
    longest = max(names)
    for i, ch in enumerate(shortest):
        if longest[i] != ch:
            return shortest[:i]
    return shortest


def get_entry_data(mibig_id: str) -> EntryData:
    filename = os.path.join(BASE_DIR, "data", f"{mibig_id}.json")
    with open(filename, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    
    publications: List[str] = data["cluster"].get("publications", [])
    compounds = [compound["compound"] for compound in data["cluster"]["compounds"]]
    name = data["cluster"]["organism_name"]
    return EntryData(publications, compounds, name)


def create_fastas(id_map: IDMap) -> Dict[str, str]:
    id_fasta: Dict[str, str] = {}
    for mibig_id, (record_id, start, end) in id_map.items():
        gbk_file = os.path.join(BASE_DIR, "genbanks", f"{record_id}.gbk")
        record = SeqIO.read(gbk_file, "genbank")
        if end == -1:
            end = len(record)
            # original start is 1-indexed
        start -= 1
        seq = record[start:end].seq
        id_fasta[mibig_id] = f">{mibig_id}|{record_id}|{start}-{end}\n{seq}"
    return id_fasta


def extract_ids(files: List[str]) -> IDMap:
    id_map: IDMap = {}
    for json_file in files:
        with open(json_file, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        if data["cluster"]["status"] == "retired":
            continue
        id_map[data["cluster"]["mibig_accession"]] = (
            data["cluster"]["loci"]["accession"],
            data["cluster"]["loci"].get("start_coord", 1),
            data["cluster"]["loci"].get("end_coord", -1),
        )

    return id_map


if __name__ == "__main__":
    main()