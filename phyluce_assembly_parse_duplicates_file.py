#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2020 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.
This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.
Created on 30 March 2020 09:44 CDT (-0500)
"""


import os
import re
import sqlite3
import argparse
import configparser
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from phyluce.helpers import (
    FullPaths,
    is_dir,
    is_file,
    get_contig_header_string
)
from phyluce.pth import get_user_param
from phyluce.log import setup_logging

import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Parse the duplicates file and output duplicate FASTA records"""
        )
    parser.add_argument(
        "--contigs",
        required=True,
        action=FullPaths,
        type=is_dir,
        help="The directory containing the assembled contigs you searched for UCE loci.",
    )
    parser.add_argument(
        '--duplicates-file',
        required=True,
        type=is_file,
        action=FullPaths,
        help='The duplicates file created by match_contigs_to_probes'
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        help="The path to the output FASTA file you want to create.",
    )
    parser.add_argument(
        "--exclude-cnt",
        type=int,
        default=0,
        help="A number of duplicate copies, above which the locus will be dropped.",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use.""",
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs.""",
    )
    return parser.parse_args()


def find_file(contigs, name):
    extensions = [
        ".fa",
        ".fasta",
        ".contigs.fasta",
        ".contigs.fa",
        ".gz",
        ".fasta.gz",
        ".fa.gz",
    ]
    for ext in extensions:
        reads1 = os.path.join(contigs, name) + ext
        reads2 = os.path.join(contigs, name.replace("-", "_")) + ext
        for reads in [reads1, reads2]:
            if os.path.isfile(reads):
                break
            elif os.path.isfile(reads.lower()):
                reads = reads.lower()
                break
            else:
                reads = None
        if reads is not None:
            break
    if reads is None:
        raise ValueError(
            "Cannot find the a fasta file for {} with any of the extensions ({}) ".format(
                name, ", ".join(extensions)
            )
        )
    return reads


def get_contig_name(header):
    """parse the contig name from the header of either velvet/trinity assembled contigs"""
    contig_header_string = get_contig_header_string()
    match = re.search("^({}).*".format(contig_header_string), header, flags=re.I)
    return match.groups()[0]


def replace_and_remove_bases(regex, seq, count):
    new_seq_string = str(seq.seq)
    if regex.search(new_seq_string):
        new_seq_string = re.sub(regex, "", new_seq_string)
        # print "\tReplaced < 20 ambiguous bases in {0}".format(seq.id)
        count += 1
    new_seq_string = re.sub("^[acgtn]+", "", new_seq_string)
    new_seq_string = re.sub("[acgtn]+$", "", new_seq_string)
    new_seq = Seq(new_seq_string)
    new_seq_record = SeqRecord(new_seq, id=seq.id, name="", description="")
    return new_seq_record, count


def main():
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    # parse the config file - allowing no values (e.g. no ":" in config file)
    config = configparser.RawConfigParser(allow_no_value=True)
    config.optionxform = str
    config.read(args.duplicates_file)
    duplicate_dict = defaultdict(dict)
    duplicate_reverse_dict = defaultdict(dict)
    # build a lookup table of which duplicate contigs go w/ which locus and taxon
    for section in config.sections():
        taxon_name = section.split(' - ')[0]
        if "probes hitting multiple contigs" in section:
            for entry in config.items(section):
                locus = entry[0]
                contigs = [i.strip() for i in entry[1].split(',')]
                for contig in contigs:
                    duplicate_dict[taxon_name][contig] = locus
                # also build a reverse list of which contigs go with which locus
                duplicate_reverse_dict[taxon_name][locus] = contigs
    # get a list of loci to drop, if duplication level is too high
    loci_to_exclude = defaultdict(set)
    if args.exclude_cnt is not 0 and args.exclude_cnt >= 1:
        for taxon in duplicate_reverse_dict.keys():
            for locus in duplicate_reverse_dict[taxon]:
                if len(duplicate_reverse_dict[taxon][locus]) > args.exclude_cnt:
                    loci_to_exclude[taxon].add(locus)
    # set up matching regex to remove Ns and ns
    regex = re.compile("[N,n]{1,21}")
    with open(args.output, "w") as uce_fasta_out:
        for organism in duplicate_dict.keys():
            # keep track of the number of contigs we associate w/ a particular locus
            duplicate_reverse_cnt_dict = defaultdict(list)
            # keep track of those loci we extract and break out of iterating over fasta file
            # if we've already got everything.
            to_get = set(duplicate_dict[organism].keys())
            gotten = set()
            text = "Getting DUPLICATE UCE loci for {0}".format(organism)
            log.info(text.center(65, "-"))
            log.info("Parsing and renaming contigs for {}".format(organism))
            name = organism.replace("_", "-")
            reads = find_file(args.contigs, name)
            count = 0
            for seq in SeqIO.parse(open(reads, "rU"), "fasta"):
                node_name = get_contig_name(seq.id)
                if node_name in duplicate_dict[organism].keys():
                    if duplicate_dict[organism][node_name] not in loci_to_exclude[organism]:
                        # add a copy to reverse dict to we can keep track of count
                        locus_name = duplicate_dict[organism][node_name]
                        duplicate_reverse_cnt_dict[locus_name].append(node_name)
                        #pdb.set_trace()
                        seq.id = "{0}_{1}_DUPE{2} |{0}".format(
                            locus_name, organism, len(duplicate_reverse_cnt_dict[locus_name])
                        )
                        seq.name = ""
                        seq.description = ""
                        # deal with strandedness because aligners sometimes dont, which
                        # is annoying
                        #if node_dict[name][1] == "-":
                        #    seq.seq = seq.seq.reverse_complement()
                        # Replace any occurrences of <21 Ns in a given sequence with
                        # blanks.  These should gap out during alignment. Also, replace
                        # leading/trailing lowercase bases from velvet assemblies.
                        # Lowercase bases indicate low coverage, and these
                        # have been problematic in downstream alignments).
                        seq, count = replace_and_remove_bases(regex, seq, count)
                        uce_fasta_out.write(seq.format("fasta"))
                    gotten.add(node_name)
                    if to_get == gotten:
                        log.info("Have reached end of list")
                        break
                    else:
                        pass
                else:
                    pass



if __name__ == "__main__":
    main()
