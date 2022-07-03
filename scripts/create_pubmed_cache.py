#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
import json
import os
import sys
from typing import List, Set, Dict

from eutils import Client

DATA_DIR = os.path.join(os.path.abspath(
    os.path.dirname(os.path.dirname(__file__))), "data")


def main():
    parser = ArgumentParser()
    parser.add_argument("cache", type=FileType(
        'w', encoding="utf-8"), help="Cache file to store pubmed cache into")
    parser.add_argument("--data-dir", type=str, default=DATA_DIR,
                        help="Directory containing the json files")
    args = parser.parse_args()

    json_files: List[str] = []
    for name in os.listdir(args.data_dir):
        if name.startswith("BGC") and name.endswith(".json"):
            json_files.append(os.path.join(args.data_dir, name))

    pmids = set()

    for json_file in json_files:
        with open(json_file) as handle:
            data = json.load(handle)
        publications = data["cluster"]["publications"]
        for publication in publications:
            if not publication.startswith("pubmed:"):
                continue
            pmids.add(publication[7:])

    collected = lookup(pmids)
    print("collected", len(collected), "entries")

    to_json = {}

    for key, val in collected.items():
        to_json[key] = val.to_json()

    json.dump(to_json, args.cache, indent=4, separators=(
        ',', ': '), sort_keys=True, ensure_ascii=False)


class PubmedEntry:
    def __init__(self, title: str, authors: List[str], year: str, journal: str, pmid: str) -> None:
        self.title = title
        self.authors = authors
        self.year = year
        self.journal = journal
        self.pmid = pmid

    @property
    def info(self) -> str:
        return f"{self.authors[0]} et al., {self.jrnl} ({self.year}) PMID:{self.pmid}"

    def to_json(self):
        return {
            "title": self.title,
            "authors": self.authors,
            "year": self.year,
            "journal": self.journal,
            "pmid": self.pmid,
        }


STEP_SIZE = 200


def lookup(pmids: Set[str]) -> Dict[str, PubmedEntry]:
    collected: Dict[str, PubmedEntry] = {}
    sorted_ids = sorted(list(pmids))
    client = Client(api_key=os.environ.get("NCBI_API_KEY", None))
    for i in range(0, len(sorted_ids), STEP_SIZE):
        articles = client.efetch(db="pubmed", id=sorted_ids[i:i+STEP_SIZE])
        for article in articles:
            collected[article.pmid] = PubmedEntry(
                article.title, article.authors, article.year, article.jrnl, article.pmid)

    return collected


if __name__ == "__main__":
    main()
