#!/usr/bin/env python3

import json
import os
from argparse import ArgumentParser

PROTEINOGENIC = {"alanine",
                 "cysteine",
                 "aspartate",
                 "glutamate"
                 "aspartic acid",
                 "glutamic acid",
                 "phenylalanine",
                 "glycine",
                 "histidine",
                 "isoleucine",
                 "lysine",
                 "leucine",
                 "methionine",
                 "asparagine",
                 "proline",
                 "glutamine",
                 "arginine",
                 "serine",
                 "threonine",
                 "valine",
                 "tryptophan",
                 "tyrosine"}

def main():
    parser = ArgumentParser("Get MIBiG database statistics")
    parser.add_argument("-old", "--old_jsons", required=True, type=str, help="Directory of previous MIBiG jsons")
    parser.add_argument("-new", "--new_jsons", required=True, type=str, help="Directory of new MIBiG jsons")
    parser.add_argument("-vold", "--version_old", default=3.0, type=float, help="Previous MIBiG version")
    parser.add_argument("-vnew", "--version_new", default=3.1, type=float, help="New MIBiG version")
    parser.add_argument("-out", "--output_dir", required=True, type=str, help="Output directory")
    args = parser.parse_args()
    run(args.old_jsons, args.new_jsons, args.version_old, args.version_new, args.output_dir)


def run(old_json_dir, new_json_dir, old_version, new_version, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    general_path = os.path.join(output_dir, "general_stats.txt")
    substrate_path_old = os.path.join(output_dir, f"substrate_stats_{old_version}.txt")
    substrate_path_new = os.path.join(output_dir, f"substrate_stats_{new_version}.txt")
    bioactivity_path_old = os.path.join(output_dir, f"bioactivity_stats_{old_version}.txt")
    bioactivity_path_new = os.path.join(output_dir, f"bioactivity_stats_{new_version}.txt")
    target_path_old = os.path.join(output_dir, f"molecular_target_stats_{old_version}.txt")
    target_path_new = os.path.join(output_dir, f"molecular_target_stats_{new_version}.txt")
    biosynthetic_classes_old = os.path.join(output_dir, f"biosyn_classes_stats_{old_version}.txt")
    biosynthetic_classes_new = os.path.join(output_dir, f"biosyn_classes_stats_{new_version}.txt")
    biosynthetic_classes_compared = os.path.join(output_dir, f"biosyn_classes_stats_compared.txt")
    bioactivity_path = os.path.join(output_dir, f"bioactivity_stats.txt")
    table = os.path.join(output_dir, f"table.txt")

    old_stats = get_stats(old_json_dir, old_version)
    new_stats = get_stats(new_json_dir, new_version)

    new_stats.removed_entries = old_stats.bgcs.difference(new_stats.bgcs)
    print(new_stats.removed_entries)
    new_stats.new_entries = new_stats.bgcs.difference(old_stats.bgcs)

    new_stats.analyse_removed_bgcs()
    new_stats.analyse_added_bgcs()

    old_stats.write_substrate_data(substrate_path_old)
    new_stats.write_substrate_data(substrate_path_new)

    old_stats.write_bioactivity_data(bioactivity_path_old)
    new_stats.write_bioactivity_data(bioactivity_path_new)

    old_stats.write_target_data(target_path_old)
    new_stats.write_target_data(target_path_new)

    old_stats.write_biosynthetic_classes(biosynthetic_classes_old)
    new_stats.write_biosynthetic_classes(biosynthetic_classes_new)
    old_stats.write_bioactivity_data_compared(new_stats, bioactivity_path)
    old_stats.write_class_data_compared(new_stats, biosynthetic_classes_compared)

    write_general_stats(old_stats, new_stats, general_path)
    write_table(new_stats, table)


def write_general_stats(stats_old, stats_new, out_file):
    with open(out_file, 'w') as out:
        out.write(f"Category\tMIBiG {stats_old.version}\tMIBiG {stats_new.version}\n")
        out.write(f"Total entries\t{stats_old.entries}\t{stats_new.entries}\n")
        out.write(f"Removed entries\t\t{len(stats_new.removed_entries)}\n")
        out.write(f"New entries\t\t{len(stats_new.new_entries)}\n")
        out.write(f"Non-minimal entries\t{stats_old.non_minimal_entries}\t{stats_new.non_minimal_entries}\n")
        out.write(f"Entries with compounds\t{stats_old.entries_with_compounds}\t{stats_new.entries_with_compounds}\n")
        out.write(f"Entries with compound structures\t{stats_old.entries_with_compound_structures}\t{stats_new.entries_with_compound_structures}\n")
        out.write(f"Entries with substrates\t{stats_old.entries_with_substrates}\t{stats_new.entries_with_substrates}\n")
        out.write(f"Entries with bioactivity\t{stats_old.entries_with_bioactivity}\t{stats_new.entries_with_bioactivity}\n")
        out.write(f"Entries with molecular targets\t{stats_old.entries_with_targets}\t{stats_new.entries_with_targets}\n")
        out.write(f"Total compounds\t{stats_old.total_compounds}\t{stats_new.total_compounds}\n")
        out.write(f"Compounds with structures\t{stats_old.compounds_with_structures}\t{stats_new.compounds_with_structures}\n")
        out.write(f"Total bioactivities\t{stats_old.total_bioactivities}\t{stats_new.total_bioactivities}\n")
        out.write(f"Compounds with bioactivities\t{stats_old.compounds_with_bioactivities}\t{stats_new.compounds_with_bioactivities}\n")
        out.write(f"Total molecular targets\t{stats_old.total_targets}\t{stats_new.total_targets}\n")
        out.write(f"Compounds with molecular targets\t{stats_old.compounds_with_molecular_targets}\t{stats_new.compounds_with_molecular_targets}\n")
        out.write(f"Compounds cross-linked to NP Atlas\t{stats_old.compounds_crosslinked_to_npatlas}\t{stats_new.compounds_crosslinked_to_npatlas}\n")
        out.write(f"Total substrates\t{stats_old.total_substrates}\t{stats_new.total_substrates}\n")
        out.write(f"Domains with substrates\t{stats_old.domains_with_substrates}\t{stats_new.domains_with_substrates}\n")
        out.write(f"Proteinogenic substrate annotations\t{stats_old.proteinogenic_substrates}\t{stats_new.proteinogenic_substrates}\n")
        out.write(f"Nonproteinogenic substrate annotations\t{stats_old.nonproteinogenic_substrates}\t{stats_new.nonproteinogenic_substrates}\n")
        out.write(f"Modular entries\t{stats_old.modular_count}\t{stats_new.modular_count}\n")


def write_table(new_stats, out_file):
    with open(out_file, 'w') as out:
        bgcs = sorted(new_stats.bgcs.union(new_stats.removed_entries))
        out.write("MIBiG accession\tNew\tRemoved\tMinimal\tHas substrate annotations\tHas compound structure\tHas bioactivity information\tUpdated substrates\tUpdated compound structure\tUpdated bioactivity information\n")
        for bgc in bgcs:
            is_new = False
            if bgc in new_stats.new_entries:
                is_new = True

            is_removed = False
            if bgc in new_stats.removed_entries:
                is_removed = True

            if is_new or is_removed:
                if bgc in new_stats.bgc_to_updates:
                    for category in new_stats.bgc_to_updates[bgc]:
                        new_stats.bgc_to_updates[bgc][category] = False
                else:
                    new_stats.bgc_to_contents[bgc] = {"bioactivity": False,
                                                      "molecular_targets": False,
                                                      "substrates": False,
                                                      "compound_structure": False}
                    new_stats.bgc_to_updates[bgc] = {"bioactivity": False,
                                                     "substrates": False,
                                                     "compound_structure": False}

            if is_removed:
                if bgc in new_stats.bgc_to_contents:
                    for category in new_stats.bgc_to_contents[bgc]:
                        new_stats.bgc_to_contents[bgc][category] = False

            out.write(f"{bgc}\t{is_new}\t{is_removed}\t{new_stats.bgc_to_minimal[bgc]}\t{new_stats.bgc_to_contents[bgc]['substrates']}\t{new_stats.bgc_to_contents[bgc]['compound_structure']}\t{new_stats.bgc_to_contents[bgc]['bioactivity']}\t{new_stats.bgc_to_updates[bgc]['substrates']}\t{new_stats.bgc_to_updates[bgc]['compound_structure']}\t{new_stats.bgc_to_updates[bgc]['bioactivity']}\n")


def get_stats(json_dir, mibig_version=3.0):
    stats = MIBiGStats(mibig_version)

    for json_file in os.listdir(json_dir):

        if json_file.endswith('.json'):

            json_path = os.path.join(json_dir, json_file)
            stats.process_json(json_path)

    return stats


class MIBiGStats:
    def __init__(self, version):
        self.version = version

        self.non_minimal_entries = 0
        self.entries = 0
        self.new_entries = set()
        self.removed_entries = set()

        self.entries_with_substrates = 0
        self.domains_with_substrates = 0
        self.nonproteinogenic_substrates = 0
        self.proteinogenic_substrates = 0
        self.total_substrates = 0

        self.bgc_to_contents = {}
        self.bgc_to_updates = {}
        self.bgc_to_minimal = {}
        self.bgc_to_class = {}

        self.entries_with_compounds = 0
        self.entries_with_compound_structures = 0
        self.total_compounds = 0

        self.compounds_with_structures = 0
        self.compounds_crosslinked_to_npatlas = 0

        self.entries_with_bioactivity = 0
        self.entries_with_targets = 0
        self.compounds_with_bioactivities = 0
        self.compounds_with_molecular_targets = 0
        self.total_bioactivities = 0
        self.total_targets = 0

        self.biosynthetic_class_to_entries = {}
        self.biosynthetic_class_to_removed_count = {}
        self.biosynthetic_class_to_added_count = {}
        self.biosynthetic_class_to_count = {}
        self.modular_count = 0

        self.bioactivity_to_count = {}
        self.target_to_count = {}
        self.substrate_to_count = {}
        self.bgcs = set()

    def process_json(self, json_path):
        with open(json_path, 'r') as json_data:

            bgc_data = json.load(json_data)
            self.bgc_to_minimal[bgc_data["cluster"]["mibig_accession"]] = True
            if self.version >= 3.0 and bgc_data['cluster']['status'] == "active":
                biosynthetic_classes = bgc_data["cluster"]['biosyn_class']
                biosynthetic_classes.sort()
                self.bgc_to_class[bgc_data["cluster"]["mibig_accession"]] = '|'.join(biosynthetic_classes)
            elif self.version < 3.0:
                biosynthetic_classes = bgc_data["cluster"]['biosyn_class']
                biosynthetic_classes.sort()
                self.bgc_to_class[bgc_data["cluster"]["mibig_accession"]] = '|'.join(biosynthetic_classes)

            if self.version >= 3.0 and bgc_data['cluster']['status'] != "active":

                if not bgc_data["cluster"]["minimal"]:
                    self.bgc_to_minimal[bgc_data["cluster"]["mibig_accession"]] = False
                return

            self.bgcs.add(bgc_data['cluster']['mibig_accession'])
            self.bgc_to_contents[bgc_data['cluster']['mibig_accession']] = {"bioactivity": False,
                                                                            "molecular_targets": False,
                                                                            "substrates": False,
                                                                            "compound_structure": False}
            self.bgc_to_updates[bgc_data['cluster']['mibig_accession']] = {"bioactivity": False,
                                                                           "substrates": False,
                                                                           "compound_structure": False}
            for changelog in bgc_data['changelog']:
                if float(changelog["version"]) == float(self.version):
                    for comment in changelog["comments"]:
                        if "Added NRP substrate specificities" in comment or "Updated NRP substrate specificities" in comment:
                            self.bgc_to_updates[bgc_data['cluster']['mibig_accession']]["substrates"] = True
                        if "Updated bioactivity data" in comment:
                            self.bgc_to_updates[bgc_data['cluster']['mibig_accession']]["bioactivity"] = True
                        if 'compound' in comment or 'structure' in comment or 'SMILES' in comment or 'smiles' in comment or 'molecule' in comment:
                            self.bgc_to_updates[bgc_data['cluster']['mibig_accession']]["compound_structure"] = True

            self.entries += 1

            if not bgc_data["cluster"]["minimal"]:
                self.non_minimal_entries += 1
                self.bgc_to_minimal[bgc_data["cluster"]["mibig_accession"]] = False

            is_polyketide = False
            counted = False

            for biosynthetic_class in sorted(biosynthetic_classes):
                if biosynthetic_class == 'NRP':
                    counted = True
                    self.modular_count += 1
                if biosynthetic_class == 'Polyketide':
                    is_polyketide = True
                if f"{biosynthetic_class}_total" not in self.biosynthetic_class_to_entries:
                    self.biosynthetic_class_to_entries[f"{biosynthetic_class}_total"] = set()
                    self.biosynthetic_class_to_count[f"{biosynthetic_class}_total"] = 0
                self.biosynthetic_class_to_entries[f"{biosynthetic_class}_total"].add(bgc_data['cluster']['mibig_accession'])
                self.biosynthetic_class_to_count[f"{biosynthetic_class}_total"] += 1
            biosynthetic_profile = '|'.join(biosynthetic_classes)
            if biosynthetic_profile not in self.biosynthetic_class_to_entries:
                self.biosynthetic_class_to_entries[biosynthetic_profile] = set()
                self.biosynthetic_class_to_count[biosynthetic_profile] = 0
            self.biosynthetic_class_to_entries[biosynthetic_profile].add(bgc_data['cluster']['mibig_accession'])
            self.biosynthetic_class_to_count[biosynthetic_profile] += 1

            entry_has_structure = False
            entry_has_bioactivity = False
            entry_has_target = False

            if is_polyketide:
                if 'polyketide' in bgc_data["cluster"] and "subclasses" in bgc_data["cluster"]["polyketide"] and "Modular type I" in bgc_data["cluster"]["polyketide"]["subclasses"]:
                    if not counted:
                        self.modular_count += 1
                        counted = True
                if 'polyketide' in bgc_data["cluster"] and "synthases" in bgc_data["cluster"]["polyketide"]:

                    for synthase in bgc_data["cluster"]["polyketide"]["synthases"]:
                        if "subclass" in synthase:
                            if "Modular type I" in synthase["subclass"]:
                                if not counted:
                                    self.modular_count += 1

            if 'compounds' in bgc_data["cluster"] and bgc_data["cluster"]["compounds"]:
                self.entries_with_compounds += 1
                compounds = bgc_data["cluster"]["compounds"]

                for compound in compounds:
                    self.total_compounds += 1
                    if "chem_struct" in compound and compound["chem_struct"]:
                        entry_has_structure = True
                        self.compounds_with_structures += 1

                    if "database_id" in compound:
                        for database_id in compound["database_id"]:
                            if database_id.startswith("npatlas"):
                                self.compounds_crosslinked_to_npatlas += 1

                    if "chem_acts" in compound and compound["chem_acts"]:
                        entry_has_bioactivity = True
                        self.compounds_with_bioactivities += 1
                        for bioactivity in compound["chem_acts"]:
                            activity = None
                            if self.version <= 3.0:
                                print(bioactivity)
                                activity = bioactivity.strip().lower()

                            elif self.version >= 3.1:
                                activity = bioactivity["activity"].strip().lower()

                            assert activity
                            if activity not in self.bioactivity_to_count:
                                self.bioactivity_to_count[activity] = 0
                            self.bioactivity_to_count[activity] += 1
                            self.total_bioactivities += 1
                    compound_has_target = False

                    if "chem_targets" in compound and compound["chem_targets"]:
                        targets = compound["chem_targets"]
                        for target in targets:
                            if 'target' in target and target['target']:
                                compound_has_target = True
                                entry_has_target = True
                                self.total_targets += 1
                                if target['target'] not in self.target_to_count:
                                    self.target_to_count[target['target']] = 0
                                self.target_to_count[target['target']] += 1

                    if compound_has_target:
                        self.compounds_with_molecular_targets += 1

            if entry_has_structure:
                self.entries_with_compound_structures += 1

            if entry_has_bioactivity:
                self.entries_with_bioactivity += 1

            if entry_has_target:
                self.entries_with_targets += 1

            entry_has_substrates = False

            if self.version >= 3.0:

                if "genes" in bgc_data["cluster"] and bgc_data["cluster"]["genes"]:
                    genes = bgc_data["cluster"]["genes"]
                    for gene in genes:
                        if "domains" in gene and gene["domains"]:
                            domains = gene["domains"]
                            for domain in domains:
                                if "name" in domain and domain["name"] == 'AMP-binding':
                                    if "substrates" in domain and domain["substrates"]:
                                        entry_has_substrates = True
                                        substrates = domain["substrates"]
                                        self.domains_with_substrates += 1
                                        for substrate in substrates:
                                            self.total_substrates += 1
                                            if substrate["name"].title() not in self.substrate_to_count:
                                                self.substrate_to_count[substrate["name"].title()] = 0
                                            self.substrate_to_count[substrate["name"].title()] += 1

                                            if substrate["name"].lower() in PROTEINOGENIC:
                                                self.proteinogenic_substrates += 1
                                            else:
                                                self.nonproteinogenic_substrates += 1
                if "nrp" in bgc_data["cluster"] and "nrps_genes" in bgc_data["cluster"]["nrp"]:
                    genes = bgc_data["cluster"]["nrp"]["nrps_genes"]
                    for gene in genes:
                        if 'modules' in gene:
                            for module in gene['modules']:
                                if 'a_substr_spec' in module and "substrates" in module['a_substr_spec']:
                                    if module['a_substr_spec']["substrates"]:
                                        substrates = module['a_substr_spec']["substrates"]
                                        entry_has_substrates = True
                                        self.domains_with_substrates += 1
                                        for substrate in substrates:
                                            if substrate["proteinogenic"]:
                                                self.proteinogenic_substrates += 1
                                            else:
                                                self.nonproteinogenic_substrates += 1
                                            self.total_substrates += 1
                                            if substrate["name"].title() not in self.substrate_to_count:
                                                self.substrate_to_count[substrate["name"].title()] = 0
                                            self.substrate_to_count[substrate["name"].title()] += 1
            elif self.version == 2.0:
                if "nrp" in bgc_data["cluster"] and "nrps_genes" in bgc_data["cluster"]["nrp"]:
                    genes = bgc_data["cluster"]["nrp"]["nrps_genes"]
                    for gene in genes:
                        if 'modules' in gene:
                            for module in gene['modules']:
                                domain_has_substrate = False
                                if 'a_substr_spec' in module:
                                    if "nonproteinogenic" in module['a_substr_spec'] and module['a_substr_spec']["nonproteinogenic"]:
                                        for substrate in module['a_substr_spec']["nonproteinogenic"]:
                                            domain_has_substrate = True
                                            entry_has_substrates = True
                                            self.nonproteinogenic_substrates += 1
                                            self.total_substrates += 1
                                            if substrate not in self.substrate_to_count:
                                                self.substrate_to_count[substrate] = 0
                                            self.substrate_to_count[substrate] += 1
                                    if "proteinogenic" in module['a_substr_spec'] and module['a_substr_spec']["proteinogenic"]:
                                        for substrate in module['a_substr_spec']["proteinogenic"]:
                                            domain_has_substrate = True
                                            entry_has_substrates = True
                                            self.proteinogenic_substrates += 1
                                            self.total_substrates += 1
                                            if substrate.title() not in self.substrate_to_count:
                                                self.substrate_to_count[substrate.title()] = 0
                                            self.substrate_to_count[substrate.title()] += 1
                                if domain_has_substrate:
                                    self.domains_with_substrates += 1

            if entry_has_substrates:
                self.entries_with_substrates += 1
                self.bgc_to_contents[bgc_data['cluster']['mibig_accession']]["substrates"] = True
            if entry_has_structure:
                self.bgc_to_contents[bgc_data['cluster']['mibig_accession']]["compound_structure"] = True
            if entry_has_bioactivity:
                self.bgc_to_contents[bgc_data['cluster']['mibig_accession']]["bioactivity"] = True
            if entry_has_target:
                self.bgc_to_contents[bgc_data['cluster']['mibig_accession']]["molecular_targets"] = True

    def analyse_added_bgcs(self):
        for biosyn_class, entries in self.biosynthetic_class_to_entries.items():
            for entry in entries:
                if entry in self.new_entries:
                    if biosyn_class not in self.biosynthetic_class_to_added_count:
                        self.biosynthetic_class_to_added_count[biosyn_class] = 0
                    self.biosynthetic_class_to_added_count[biosyn_class] += 1

    def analyse_removed_bgcs(self):
        for biosyn_class, entries in self.biosynthetic_class_to_entries.items():
            for entry in entries:
                if entry in self.removed_entries:
                    if biosyn_class not in self.biosynthetic_class_to_removed_count:
                        self.biosynthetic_class_to_removed_count[biosyn_class] = 0
                    self.biosynthetic_class_to_removed_count[biosyn_class] += 1

    def write_biosynthetic_classes(self, out_path):
        with open(out_path, 'w') as out:
            out.write('Biosynthetic class\tCount\n')
            for biosyn_class, count in self.biosynthetic_class_to_count.items():
                out.write(f"{biosyn_class}\t{count}\n")

    def write_substrate_data(self, out_path):
        with open(out_path, 'w') as out:
            out.write('Substrate\tCount\n')
            for substrate, count in self.substrate_to_count.items():
                out.write(f"{substrate}\t{count}\n")

    def write_target_data(self, out_path):
        with open(out_path, 'w') as out:
            out.write('Molecular target\tCount\n')
            for target, count in self.target_to_count.items():
                out.write(f"{target}\t{count}\n")

    def write_bioactivity_data(self, out_path):
        with open(out_path, 'w') as out:
            out.write('Bioactivity\tCount\n')
            for bioactivity, count in self.bioactivity_to_count.items():
                out.write(f"{bioactivity}\t{count}\n")

    def write_class_data_compared(self, other, out_path):
        old_bgcs = set()
        new_bgcs = set()
        retired_bgcs = set()
        changed_bgcs = set()
        class_to_figures = {}

        for bgc, biosyn_class in other.bgc_to_class.items():
            if biosyn_class not in class_to_figures:
                class_to_figures[biosyn_class] = {"old": 0,
                                                  "new": 0,
                                                  "retired": 0,
                                                  "changed": 0}
                for subclass in biosyn_class.split('|'):
                    if f"{subclass}_total" not in class_to_figures:
                        class_to_figures[f"{subclass}_total"] = {"old": 0,
                                                                 "new": 0,
                                                                 "retired": 0,
                                                                 "changed": 0}

            if bgc in self.bgc_to_class:
                if biosyn_class == self.bgc_to_class[bgc]:
                    old_bgcs.add(bgc)
                    class_to_figures[biosyn_class]["old"] += 1
                    for subclass in biosyn_class.split('|'):
                        class_to_figures[f"{subclass}_total"]["old"] += 1
                else:
                    changed_bgcs.add(bgc)
                    class_to_figures[biosyn_class]["changed"] += 1
                    for subclass in biosyn_class.split('|'):
                        class_to_figures[f"{subclass}_total"]["changed"] += 1

            else:
                new_bgcs.add(bgc)
                class_to_figures[biosyn_class]["new"] += 1
                for subclass in biosyn_class.split('|'):
                    class_to_figures[f"{subclass}_total"]["new"] += 1

        for bgc, biosyn_class in self.bgc_to_class.items():
            if biosyn_class not in class_to_figures:
                class_to_figures[biosyn_class] = {"old": 0,
                                                  "new": 0,
                                                  "retired": 0,
                                                  "changed": 0}

            if bgc not in other.bgc_to_class:
                retired_bgcs.add(bgc)
                class_to_figures[biosyn_class]["retired"] += 1
                for subclass in biosyn_class.split('|'):
                    class_to_figures[f"{subclass}_total"]["retired"] += 1

        with open(out_path, 'w') as out:
            out.write(f"Biosynthetic class\tKept from {self.version}\tNew in {other.version}\tChanged in {other.version}\tRetired in {other.version}\n")
            out.write(f"All\t{len(old_bgcs)}\t{len(new_bgcs)}\t{len(changed_bgcs)}\t{len(retired_bgcs)}\n")
            for biosyn_class, figures in class_to_figures.items():
                out.write(f"{biosyn_class}\t{figures['old']}\t{figures['new']}\t{figures['changed']}\t{figures['retired']}\n")

    def write_bioactivity_data_compared(self, other, out_path):
        seen_categories = set()
        with open(out_path, 'w') as out:
            out.write(f'Bioactivity\tCount {self.version}\tCount {other.version}\n')
            for bioactivity, count in self.bioactivity_to_count.items():
                if bioactivity.lower() in other.bioactivity_to_count:
                    out.write(f"{bioactivity}\t{count}\t{other.bioactivity_to_count[bioactivity.lower()]}\n")
                else:
                    out.write(f"{bioactivity}\t{count}\t{0}\n")
                seen_categories.add(bioactivity.lower())
            for bioactivity, count in other.bioactivity_to_count.items():
                if bioactivity.lower() not in seen_categories:
                    out.write(f"{bioactivity}\t{0}\t{other.bioactivity_to_count[bioactivity.lower()]}\n")


if __name__ == "__main__":
    main()

