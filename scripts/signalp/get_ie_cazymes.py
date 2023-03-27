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
"""Add intracellular/extracellular annotations to CAZy family annotations"""

import pandas as pd

from saintBioutils.utilities.file_io.get_paths import get_file_paths
from tqdm import tqdm


SIGNALP_FILE="data/pecto_dic/signalp/pd_signalp_output"
THRESHOLD=0.95
FGP_FILE="data/pecto_dic/cazomes/pd_fam_genomes_proteins"

print("Loading CAZy family annotations into dict")
fgp_df = pd.read_table(FGP_FILE, header=None)
fgp_df.columns = ['Fam', 'Genome', 'Protein']

print("Parse signalP output, and identify intra/extracellular CAZymes")
signalp_df = pd.read_table(SIGNALP_FILE, header=1)
signalp_df = signalp_df[['# ID', 'CS Position']]

signalp_output = {
    "intra": set(),  # sets of protein accs - intracellular
    "extra": set(),  # extracellular
}
for ri in tqdm(range(len(signalp_df)), desc="Identifying intra- and extra-cellular CAZymes"):
    protein_id = signalp_df.iloc[ri]['# ID'].split(" ")[0]
    intracellular = True

    if str(signalp_df.iloc[ri]['CS Position']).strip() != 'nan':
        if float(signalp_df.iloc[ri]['CS Position'].split(" ")[-1]) >= THRESHOLD:
            intracellular = False

    if intracellular:
        signalp_output['intra'].add(protein_id)
    else:
        signalp_output['extra'].add(protein_id)

ie_data = []
for ri in tqdm(range(len(fgp_df)), desc="Adding intra-/extra-cellular annotations to families"):
    fam = fgp_df.iloc[ri]['Fam']
    genome = fgp_df.iloc[ri]['Genome']
    protein_id = fgp_df.iloc[ri]['Protein']

    if protein_id in signalp_output['intra']:
        fam = f"i_{fam}"
    else:
        fam = f"e_{fam}"
    
    ie_data.append([fam, genome, protein_id])

ie_fgp_df = pd.DataFrame(ie_data, columns=['Fam','Genome','Protein'])
ie_fgp_df.to_csv("data/pecto_dic/cazomes/pd_fam_genomes_proteins", sep="\t")
