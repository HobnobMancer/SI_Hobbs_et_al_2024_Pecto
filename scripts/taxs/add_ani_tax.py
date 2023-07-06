#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Add tax data to ANI tree"""


from tqdm import tqdm


TREE_FILE = "data/tree/ani_tree/pyani_ani_tree.new"
TAB_FILE = "data/cazomes/pecto_fam_genomes_taxs"
OUTFILE = "data/tree/ani_tree/pyani_ani_tree_taxs.new"


with open(TREE_FILE, 'r') as fh:
    tree = fh.read()

genome_dict = {}  # {genome: 'genome_tax'}

with open(TAB_FILE, 'r') as fh:
    lines = fh.read().splitlines()

for line in tqdm(lines, desc="Parsing tab delimited list of genome_tax and fam annotations"):
    genome = line.split("\t")[1]
    genome = f"{genome.split('_')[0]}_{genome.split('_')[1]}"
    genome_dict[genome] = line.split("\t")[1].replace(" ",'_')

# replace genomes in tree
for genome in tqdm(genome_dict, desc="Adding tax data to genomes"):
    tree = tree.replace(genome, genome_dict[genome])

with open(OUTFILE, 'w') as fh:
    fh.write(tree)
