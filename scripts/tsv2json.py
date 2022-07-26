#!/usr/bin/env python3

from argparse import (
    ArgumentParser,
    FileType,
)
from collections import defaultdict
import json
import os
import re
import sys
from typing import (
    Dict,
    IO,
    List,
    Optional,
    Tuple,
    Union,
)

from ratelimiter import RateLimiter
import requests
import untangle

import pubchempy as pcp


MAX_CALLS = 3
MAX_CALLS_NPATLAS = 1
SLASH_SEP = re.compile(r' / ')

StringOrInt = Union[str, int]
StringOrStringList = Union[str, List[str]]

TaxCache = Dict[str, Dict[str, StringOrInt]]
SmilesCache = Dict[str, str]
TCache = Dict[str, Union[TaxCache, SmilesCache]]
ChangeLog = List[Dict[str, StringOrStringList]]
TLocus = Dict[str, Union[str, List[str], int]]
TCompound = Dict[str, StringOrStringList]
Cluster = Dict[str, Union[StringOrStringList, TLocus, List[TCompound], bool, int]]


def main():
    parser = ArgumentParser("Parse the MIBiG 3.0 TSV dataset into JSON")
    parser.add_argument("infile", type=FileType("r", encoding="utf-8"),
                        help="TSV file to parse")
    parser.add_argument("-o", "--outdir", default="data",
                        help="Directory to create JSON files in")
    parser.add_argument("-c", "--cache", default="cache.json",
                        help="Cache for the downloaded info")
    parser.add_argument("-s", "--structures", type=str, default='',
                        help="TSV file containing chemical structure information")
    parser.add_argument("-np", "--npatlas_database", type=str, default='')

    args = parser.parse_args()

    run(args.infile, args.outdir, args.cache,
        args.structures, args.npatlas_database)


def run(infile: IO, outdir: str, cache_file: str, structure_file: str = '',
        npatlas_database: str = '') -> None:

    cache = Cache(cache_file)

    seen_ids = set()

    bgc_to_compounds: Dict[str, List[StructureEntry]] = defaultdict(list)
    npatlas_data: Optional[NPAtlasData] = None

    if npatlas_database:
        npatlas_data = NPAtlasData(npatlas_database)
        cache.load_smiles_from_npatlas(npatlas_data)

    if structure_file:
        with open(structure_file, 'r') as structure_data:
            for structure_entry in StructureReader(structure_data):
                bgc_to_compounds[structure_entry.bgc_accession].append(structure_entry)

    for entry in EntryReader(infile, cache):
        if entry.accession in seen_ids:
            print(entry.accession, "is a duplicate", file=sys.stderr)
        else:
            seen_ids.add(entry.accession)
            if entry.tmp_accession in bgc_to_compounds:
                non_matching_compounds = []
                matching_compounds = []
                for compound in entry.compounds:
                    lc_names = (c.name.lower() for c in bgc_to_compounds[entry.tmp_accession])
                    if compound.name.lower() not in lc_names:
                        non_matching_compounds.append(compound)
                    else:
                        matching_compounds.append(compound)
                if matching_compounds:
                    entry.compounds = non_matching_compounds
                    if non_matching_compounds:
                        unmatched = ', '.join([c.name for c in non_matching_compounds])
                        print(f"Review the following compounds in {entry.accession}: {unmatched}")
                else:
                    entry.compounds = []
                    skipped = ', '.join([c.name for c in non_matching_compounds])
                    print(f"Skipped the following compounds in {entry.accession}: {skipped}.")

                for structure in bgc_to_compounds[entry.tmp_accession]:
                    smiles = ''
                    if structure.npaid and npatlas_data:
                        smiles = npatlas_data.get_part(
                            'compound_smiles', structure.npaid)
                    compound = Compound(structure.name, cache, structure.npaid, structure.pid,
                                        structure.chemspider,
                                        structure.lotus, structure.chembl, smiles)
                    for original_compound in matching_compounds:
                        if original_compound.name.lower() == structure.name.lower():
                            if original_compound.synonyms:
                                compound.synonyms = original_compound.synonyms

                    for synonym in structure.synonyms:
                        if synonym.lower() not in [sym.lower() for sym in compound.synonyms]:
                            compound.synonyms.append(synonym)

                    if not compound.smiles:
                        compound.smiles = structure.smiles
                    entry.compounds.append(compound)
            else:
                for compound in entry.compounds:
                    if compound.npaid and npatlas_data:
                        smiles = npatlas_data.get_part(
                            'compound_smiles', compound.npaid)
                        compound.smiles = smiles

        name = f"{entry.accession}.json"
        if " " in name:
            name = name.replace(" ", "_")
        with open(os.path.join(outdir, name), "w", encoding="utf-8") as handle:
            json.dump(entry.to_json(), handle, indent=4, ensure_ascii=False)


class StructureEntry:
    def __init__(self, headers: List[str], line: str) -> None:
        self._headers = headers
        self._parts = line.rstrip('\n').split('\t')
        self.bgc_accession = self._get_part('bgc_id').strip()
        if self.bgc_accession.startswith("NEW0"):
            self.bgc_accession = "BGC9" + self.bgc_accession[4:]
        self.name = self._get_part('compound_name').strip()
        self.synonyms = []
        if ' = ' in self.name:
            names = self.name.split(' = ')
            self.name = names[0]
            self.synonyms = names[1:]
            self.synonyms.sort()
        self.npaid = self._get_part('npaid').strip()
        self.lotus = self._get_part('lotus').strip()
        self.pid = self._get_part('pid').strip()
        self.chemspider = self._get_part('chemspider').strip()
        self.chembl = self._get_part('chembl').strip()
        self.smiles = self._get_part('SMILES').strip()

        if self.chembl and self.chembl[:6].lower() == 'chembl':
            self.chembl = 'CHEMBL' + self.chembl[6:]
        elif self.chembl:
            self.chembl = 'CHEMBL' + self.chembl

    def _get_part(self, name: str) -> str:
        try:
            idx = self._headers.index(name)
        except ValueError:
            print(name, "not in", self._headers)
            raise

        return self._parts[idx]


class Entry:
    def __init__(self, headers: List[str], line: str, cache: "Cache") -> None:
        self._headers = headers
        self._parts = line.rstrip('\n').split('\t')
        # Google Docs doesn't add the empty comment field...
        if len(self._parts) != len(self._headers) and len(self._parts) != len(self._headers) - 1:
            acc = ""
            if len(self._parts) > 1:
                acc = self._parts[1] + ": "
            min_len = min(len(self._headers), len(self._parts))
            comp = ""
            for i in range(min_len):
                comp += f"{self._headers[i]} <> {self._parts[i]}\n"
            raise ValueError(
                f"{acc}Headers don't match parts:\n{comp}")
        self.accession = self._get_part("Final ID")
        if self.accession.startswith("NEW0"):
            raise ValueError(f"invalid accession {self.accession}")

        self.tmp_accession = self._get_part("temporary id")
        if self.tmp_accession.startswith('NEW0'):
            self.tmp_accession = 'BGC9' + self.tmp_accession[4:]

        bgc_types = [self._get_part("biosynthetic class")]
        if self._get_part("biosyn class 2"):
            bgc_types.append(self._get_part("biosyn class 2"))
        if self._get_part("biosyn class 3"):
            bgc_types.append(self._get_part("biosyn class 3"))
        self.biosynthetic_classes = bgc_types

        compound_fields = self._get_part("compound name(s)").split("/")
        npaid_fields = self._get_part("NPAID").split("/")

        self.compounds: List[Compound] = []
        has_npaid = False

        if len(compound_fields) == len(npaid_fields):
            has_npaid = True
        for i, compound in enumerate(compound_fields):
            compound = compound.strip()
            if not has_npaid:
                self.compounds.append(Compound(compound, cache))
            else:
                npaid = npaid_fields[i].strip()
                if npaid.startswith('NPA'):
                    self.compounds.append(
                        Compound(compound, cache, npaid=npaid))
                else:
                    self.compounds.append(Compound(compound, cache))

        db_acc = self._get_part("GenBank/RefSeq acc").strip()
        start = self._get_part("start coord")
        end = self._get_part("end coord")
        evidence = [self._get_part("evidence 1"), self._get_part(
            "evidence 2"), self._get_part("evidence 3")]

        try:
            self.organism_name, self.taxid, db_acc = cache.fetch_taxon_info(db_acc)
        except ValueError as err:
            raise ValueError(f"{self.accession}: {err}")

        self.locus = Locus(db_acc, start, end, evidence)

        self.publications = Publications(self._get_part("PubMed ID"), self._get_part(
            "doi"))

    def _get_part(self, name: str) -> str:
        try:
            idx = self._headers.index(name)
        except ValueError:
            print(name, "not in", self._headers)
            raise

        return self._parts[idx]

    def to_json(self) -> Dict[str, Union[ChangeLog, Cluster]]:
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
                "status": "active",
            },
        }


class Compound:
    def __init__(self, name: str, cache: "Cache", npaid: str = '', pid: str = '',
                 chemspider: str = '', lotus: str = '', chembl: str = '',
                 smiles: str = '') -> None:
        self.name = name
        self.synonyms = []
        if ' = ' in name:
            names = name.split(' = ')
            self.name = names[0]
            self.synonyms = names[1:]
            self.synonyms.sort()

        self.npaid = npaid.strip()
        self.pid = pid.strip()
        self.smiles = ''

        if smiles:
            self.smiles = smiles
        elif self.npaid:
            self.smiles = cache.get_smiles(self.npaid)
        elif self.pid:
            self.smiles = cache.get_smiles(self.pid)

        self.chemspider = chemspider
        self.lotus = lotus
        self.chembl = chembl

    def to_json(self) -> TCompound:
        json_dict: TCompound = {"compound": self.name}
        if self.npaid or self.pid or self.lotus or self.chembl or self.chemspider:
            db_ids = []
            if self.npaid:
                db_ids.append(f'npatlas:{self.npaid}')
            if self.pid:
                db_ids.append(f'pubchem:{self.pid}')
            if self.lotus:
                db_ids.append(f'lotus:{self.lotus}')
            if self.chembl:
                db_ids.append(f'chembl:{self.chembl}')
            if self.chemspider:
                db_ids.append(f'chemspider:{self.chemspider}')

            if db_ids:
                json_dict["database_id"] = db_ids

        if self.smiles:
            json_dict.update({"chem_struct": self.smiles})
        if self.synonyms:
            json_dict["chem_synonyms"] = self.synonyms

        return json_dict


class Locus:
    def __init__(self, acc: str, start: str, end: str, evidence: List[str]) -> None:
        self.completeness = "Unknown"
        self.evidence = [e for e in evidence if e and e != "Other"]
        self.acc = acc
        try:
            self.start: Optional[int] = int(start)
        except ValueError:
            self.start = None

        try:
            self.end: Optional[int] = int(end)
        except ValueError:
            self.end = None

    def to_json(self) -> TLocus:
        ret: TLocus = {
            "accession": self.acc,
            "completeness": self.completeness,
            "evidence": self.evidence,
        }
        if self.start:
            ret["start_coord"] = self.start
        if self.end:
            ret["end_coord"] = self.end

        return ret


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


class NPAtlasData:
    def __init__(self, tsv_file: str):
        self.data: Dict[str, List[str]] = {}
        with open(tsv_file, 'r') as tsv:
            reader = TsvReader(tsv)
            for line in reader:
                line = line.strip('\n')
                if line:
                    data = line.split('\t')
                    row_id = data[0]
                    self.data[row_id] = data
            self.categories = reader.headers

    def get_part(self, name: str, row_id: str) -> str:
        try:
            idx = self.categories.index(name)
        except ValueError:
            print(name, "not in", self.categories)
            raise

        try:
            row = self.data[row_id]
        except KeyError:
            print(row_id, "not in data.")
            raise

        return row[idx]


class TsvReader:
    def __init__(self, infile: IO) -> None:
        self.infile = infile
        self.count = 0

    def __iter__(self) -> "TsvReader":
        # Read the first line to get the header fields
        self.line = self.infile.readline()
        self.headers = self.line.strip().split('\t')
        return self

    def __next__(self) -> str:
        line = self.infile.readline()
        self.count += 1
        if not line:
            raise StopIteration
        return line


class StructureReader:
    def __init__(self, infile: IO) -> None:
        self._reader = TsvReader(infile)

    def __iter__(self) -> "StructureReader":
        self._reader.__iter__()
        self.headers = self._reader.headers
        return self

    def __next__(self) -> StructureEntry:
        line = next(self._reader)
        try:
            entry = StructureEntry(self.headers, line)
            return entry
        except ValueError as err:
            print("error in line", self._reader.count, err, file=sys.stderr)
            return next(self)


class EntryReader:
    def __init__(self, infile: IO, cache: "Cache") -> None:
        self._reader = TsvReader(infile)
        self.cache = cache

    def __iter__(self) -> "EntryReader":
        self._reader.__iter__()
        self.headers = self._reader.headers
        return self

    def __next__(self) -> Entry:
        line = next(self._reader)
        try:
            entry = Entry(self.headers, line, self.cache)
            return entry
        except ValueError as err:
            print("error in line", self._reader.count, err, file=sys.stderr)
            return next(self)


class Cache:
    def __init__(self, filename: str) -> None:
        self.filename = filename
        self.tax_cache: TaxCache = {}
        self.smiles_cache: SmilesCache = {}
        if os.path.exists(filename):
            with open(filename, 'r') as handle:
                cache = json.load(handle)
            self.tax_cache = cache["tax_cache"]
            self.smiles_cache = cache["smiles_cache"]

        # side-load one missing taxon entry
        self.tax_cache["ON409579"] = {
            "orgname": "Candidatus Thermopylae lasonolidus",
            "taxid": 134622,
            "accession": "ON409579.1"
        }

    def fetch_taxon_info(self, accession: str) -> Tuple[str, int, str]:
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

        if accession in self.tax_cache:
            return (str(self.tax_cache[accession]["orgname"]),
                    int(self.tax_cache[accession]["taxid"]),
                    str(self.tax_cache[accession]['accession']))

        raw = _get_from_ncbi(accession)
        data = untangle.parse(raw)
        try:
            orgname = data.TSeqSet.TSeq.TSeq_orgname.cdata
            taxid = int(data.TSeqSet.TSeq.TSeq_taxid.cdata)
            accession_versioned = data.TSeqSet.TSeq.TSeq_accver.cdata
            self.tax_cache[accession] = {}
            self.tax_cache[accession]["orgname"] = orgname
            self.tax_cache[accession]["taxid"] = taxid
            self.tax_cache[accession]["accession"] = accession_versioned

            self.save()

            return orgname, taxid, accession_versioned
        except AttributeError as err:
            raise ValueError(f"{accession} raised error {err}")

    def load_smiles_from_npatlas(self, npatlas_data: NPAtlasData) -> None:
        for key in npatlas_data.data.keys():
            self.smiles_cache[key] = npatlas_data.get_part('compound_smiles', key)
        self.save()

    def get_smiles(self, compound_id: str) -> str:
        if compound_id in self.smiles_cache:
            return self.smiles_cache[compound_id]

        if compound_id.startswith("NPA"):
            smiles = _get_smiles_npatlas(compound_id)
        else:
            smiles = _get_smiles_pubchem(compound_id)

        self.smiles_cache[compound_id] = smiles

        self.save()

        return smiles

    def save(self) -> None:
        cache = {
            "tax_cache": self.tax_cache,
            "smiles_cache": self.smiles_cache
        }
        with open(self.filename, 'w', encoding="utf-8") as handle:
            json.dump(cache, handle, ensure_ascii=False)


@RateLimiter(MAX_CALLS)
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


@RateLimiter(MAX_CALLS_NPATLAS)
def _get_smiles_npatlas(npaid: str) -> str:
    try:
        URL = f'https://www.npatlas.org/api/v1/compound/{npaid}?include=classifications'

        response = requests.get(URL)
        json_data = json.loads(response.content)
        return json_data['smiles']
    except Exception as err:
        raise ValueError(f"Could not find SMILES for NP Atlas ID {npaid}: {err}.")


@RateLimiter(MAX_CALLS)
def _get_smiles_pubchem(pubchem_id_str: str) -> str:
    try:
        pubchem_id = int(pubchem_id_str)
        c = pcp.Compound.from_cid(pubchem_id)
        if c.isomeric_smiles:
            return str(c.isomeric_smiles)
        else:
            return str(c.canonical_smiles)
    except Exception as err:
        raise ValueError(f"Could not find SMILES for pubchem ID {pubchem_id}: {err}.")


if __name__ == "__main__":
    main()
