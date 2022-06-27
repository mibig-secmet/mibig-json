#!/usr/bin/env python3

from argparse import (
    ArgumentParser,
    FileType,
)
import json
import os
import re
import sys
from typing import (
    Any,
    Dict,
    IO,
    List,
    Tuple,
    Union,
)

from ratelimiter import RateLimiter
import requests
import untangle


MAX_CALLS = 3

Cache = Dict[str, Dict[str, Union[str, int]]]


def main():
    parser = ArgumentParser("Parse the MIBiG 3.0 TSV dataset into JSON")
    parser.add_argument("infile", type=FileType(
        "r", encoding="utf-8"), help="TSV file to parse")
    parser.add_argument("-o", "--outdir", default="data",
                        help="Directory to create JSON files in")
    parser.add_argument(
        "-c", "--cache", default="ncbi_cache.json", help="Cache for the NCBI info")
    args = parser.parse_args()
    run(args.infile, args.outdir, args.cache)


def run(infile: IO, outdir: str, cache_file: str) -> None:

    cache: Cache = {}

    if os.path.exists(cache_file):
        with open(cache_file, 'r') as handle:
            cache = json.load(handle)

    seen_ids = set()

    for entry in TsvReader(infile, cache):
        if entry.accession in seen_ids:
            print(entry.accession, "is a duplicate", file=sys.stderr)
        else:
            seen_ids.add(entry.accession)
        name = f"{entry.accession}.json"
        if " " in name:
            name = name.replace(" ", "_")
        with open(os.path.join(outdir, name), "w", encoding="utf-8") as handle:
            json.dump(entry.to_json(), handle, indent=4, ensure_ascii=False)

    with open(cache_file, 'w', encoding="utf-8") as handle:
        json.dump(cache, handle, ensure_ascii=False)


class Entry:
    def __init__(self, headers: List[str], line: str, cache: Cache) -> None:
        self._headers = headers
        self._parts = line.rstrip('\n').split('\t')
        # Google Docs doesn't add the empty comment field...
        if len(self._parts) != len(self._headers) and len(self._parts) != len(self._headers) - 1:
            acc = ""
            if len(self._parts) > 1:
                acc = self._parts[1] + ": "
            l = min(len(self._headers), len(self._parts))
            comp = ""
            for i in range(l):
                comp += f"{self._headers[i]} <> {self._parts[i]}\n"
            raise ValueError(
                f"{acc}Headers don't match parts:\n{comp}")
        self.accession = self._get_part("temporary id")
        if self.accession.startswith("NEW0"):
            self.accession = "BGC9" + self.accession[4:]
        bgc_types = [self._get_part("biosynthetic class")]
        if self._get_part("biosyn class 2"):
            bgc_types.append(self._get_part("biosyn class 2"))
        if self._get_part("biosyn class 3"):
            bgc_types.append(self._get_part("biosyn class 3"))
        self.biosynthetic_classes = bgc_types
        compound_fields = self._get_part("compound name(s)").split("/")
        self.compounds = []
        for compound in compound_fields:
            self.compounds.append(Compound(compound.strip()))

        db_acc = self._get_part("GenBank/RefSeq acc").strip()
        start = self._get_part("start coord")
        end = self._get_part("end coord")
        evidence = [self._get_part("evidence 1"), self._get_part(
            "evidence 2"), self._get_part("evidence 3")]
        self.locus = Locus(db_acc, start, end, evidence)

        try:
            self.organism_name, self.taxid = fetch_taxon_info(db_acc, cache)
        except ValueError as err:
            raise ValueError(f"{self.accession}: {err}")

        self.publications = Publications(self._get_part("PubMed ID"), self._get_part(
            "doi"))

    def _get_part(self, name: str) -> str:
        try:
            idx = self._headers.index(name)
        except ValueError:
            print(name, "not in", self._headers)
            raise

        return self._parts[idx]

    def to_json(self) -> Dict[str, Union[str, Dict[str, Union[str, int]]]]:
        return {
            "changelog": [{
                "comments": ["Entry added"],
                "contributors": ["AAAAAAAAAAAAAAAAAAAAAAAA"],
                "version": "3.0",
            }],
            "cluster": {
                "biosyn_class": self.biosynthetic_classes,
                "compounds": [c.to_json() for c in self.compounds],
                "loci": self.locus.to_json(),
                "mibig_accession": self.accession,
                "minimal": True,
                "organism_name": self.organism_name,
                # TODO: Make taxid a number in the schema
                "ncbi_tax_id": f"{self.taxid}",
                "publications": self.publications.to_json(),
            },
        }


class Compound:
    def __init__(self, name: str) -> None:
        self.name = name

    def to_json(self) -> Dict[str, str]:
        return {"compound": self.name}


class Locus:
    def __init__(self, acc: str, start: str, end: str, evidence: List[str]) -> None:
        self.completeness = "Unknown"
        self.evidence = [e for e in evidence if e and e != "Other"]
        self.acc = acc
        try:
            self.start = int(start)
        except ValueError:
            self.start = None

        try:
            self.end = int(end)
        except ValueError:
            self.end = None

    def to_json(self) -> Dict[str, Union[str, int, List[str]]]:
        ret = {
            "accession": self.acc,
            "completeness": self.completeness,
            "evidence": self.evidence,
        }
        if self.start:
            ret["start_coord"] = self.start
        if self.end:
            ret["end_coord"] = self.end

        return ret


SLASH_SEP = re.compile(r' / ')


class Publications:
    def __init__(self, pmid: str, doi: str) -> None:
        self.pmid = pmid
        self.doi = doi

    def to_json(self) -> List[str]:
        ret: List[str] = []
        if self.pmid:
            pmids: List[str] = []
            pmids = [p.strip() for p in re.split(SLASH_SEP, self.pmid)]
            if not pmids:
                pmids = [f"pubmed:{self.pmid}"]
            for pmid in pmids:
                if not pmid.isnumeric():
                    continue
                pmid = pmid.split(":")[-1].strip()
                ret.append(f"pubmed:{pmid}")
        if self.doi:
            dois: List[str] = []
            dois = [d.strip() for d in re.split(SLASH_SEP, self.doi)]
            for doi in dois:
                doi = self._clean_doi(doi)
                ret.append(f"doi:{doi}")

        return ret

    def _clean_doi(self, doi: str) -> str:
        if doi.startswith("https://doi.org/"):
            return doi[16:]
        if doi.startswith("http://doi.org/"):
            return doi[15:]
        if doi.startswith("doi.org/"):
            return doi[8:]
        if doi.startswith("DOI:") or doi.startswith("doi:") or doi.startswith("doi/"):
            return doi[4:].strip()
        if doi.startswith("https://doi-org.proxy-ub.rug.nl/"):
            return doi[32:]
        return doi


class TsvReader:
    def __init__(self, infile: IO, cache: Cache) -> None:
        self.infile = infile
        self.cache = cache
        self.count = 0

    def __iter__(self) -> "TsvReader":
        # Read the first line to get the header fields
        self.line = self.infile.readline()
        self.headers = self.line.strip().split('\t')
        return self

    def __next__(self) -> Entry:
        self.line = self.infile.readline()
        self.count += 1
        if not self.line:
            raise StopIteration
        try:
            entry = Entry(self.headers, self.line, self.cache)
            return entry
        except ValueError as err:
            print("error in line", self.count, err, file=sys.stderr)
            return next(self)


def fetch_taxon_info(accession: str, cache: Cache) -> Tuple[str, int]:
    if not accession:
        raise ValueError("NCBI accession is empty")
    if " " in accession:
        raise ValueError(f"{accession} contains a space")
    if "," in accession:
        raise ValueError(f"{accession} contains a comma")
    if "-" in accession:
        raise ValueError(f"{accession} contains a dash")
    if "|" in accession:
        raise ValueError(f"{accession} contains a pipe")
    if "/" in accession:
        raise ValueError(f"{accession} contains a slash")
    if accession.startswith("GCA_") or accession.startswith("GCF_"):
        raise ValueError(f"{accession} is an assembly ID")
    if len(accession) < 5:
        raise ValueError(f"{accession} is too short")
    if accession.startswith("WP_") or accession.startswith("YP_"):
        raise ValueError(f"{accession} is a protein ID")
    if accession.split(".")[0].endswith("000000"):
        raise ValueError(
            f"{accession} is a WGS record, please supply the actual contig")
    if accession.startswith("PRJ"):
        raise ValueError(
            f"{accession} is a bioproject id, please supply actual GenBank record")
    if accession.startswith("SRR"):
        raise ValueError(
            f"{accession} is a SRA record, please supply an assembled contig")

    if accession in cache:
        return cache[accession]["orgname"], cache[accession]["taxid"]

    raw = _get_from_ncbi(accession)
    data = untangle.parse(raw)
    try:
        orgname = data.TSeqSet.TSeq.TSeq_orgname.cdata
        taxid = int(data.TSeqSet.TSeq.TSeq_taxid.cdata)
        cache[accession] = {}
        cache[accession]["orgname"] = orgname
        cache[accession]["taxid"] = taxid

        return orgname, taxid
    except AttributeError as err:
        raise ValueError(f"{accession} raised error {err}")


@ RateLimiter(MAX_CALLS)
def _get_from_ncbi(accession: str) -> str:
    params = {
        "db": "nuccore",
        "rettype": "fasta",
        "retmode": "xml",
        "id": accession,
        "seq_start": "1",
        "seq_stop": "1",
    }
    res = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params)

    if res.status_code != 200:
        raise ValueError(f"{accession}: NCBI query returned {res.status_code}")

    return res.text


if __name__ == "__main__":
    main()
