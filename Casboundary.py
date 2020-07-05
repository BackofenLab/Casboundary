"""
    Casboundary
    Copyright (C) 2020 Victor Alexandre Padilha <victorpadilha@usp.br>,
                       Omer Salem Alkhnbashi <alkhanbo@informatik.uni-freiburg.de>,
                       Van Dinh Tran <dinh@informatik.uni-freiburg.de>,
                       Shiraz Ali Shah <shiraz.shah@dbac.dk>,
                       André Carlos Ponce de Leon Ferreira de Carvalho <andre@icmc.usp.br>,
                       Rolf Backofen <backofen@informatik.uni-freiburg.de>
    
    This file is part of Casboundary.
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os
import glob
import numpy as np
import warnings
warnings.simplefilter('ignore')

from pathlib import Path

# Project imports
from prodigal import run_prodigal, prodigal_fasta_to_genome_info
from hmmer import run_hmmsearch
from utils import extract_targz

CAS_HMM_TAR_PATH = 'Cas_HMM.tar.gz'
SIG_HMM_TAR_PATH = 'Sig_HMM.tar.gz'
GEN_HMM_TAR_PATH = 'Gen_HMM.tar.gz'

CAS_HMM_DIR = 'Cas_HMM'
SIG_HMM_DIR = 'Sig_HMM'
GEN_HMM_DIR = 'Gen_HMM'

GEN_HMM_NAMES = 'gen_hmm_names.txt'

def find_potential_regions(fasta_file, sequence_completeness, output_dir, hmmsearch_output_dir, n_cpus, offset=50):
    print('Running Prodigal on input Fasta file.')
    proteins_fasta_file = run_prodigal(fasta_file, sequence_completeness, output_dir)

    print('Generating Genome DataFrame')
    info_df = prodigal_fasta_to_genome_info(proteins_fasta_file, output_dir)

    print('Searching for potential signature proteins')
    sig_output_dir = hmmsearch_output_dir + '/' + SIG_HMM_DIR
    run_hmmsearch(proteins_fasta_file, SIG_HMM_DIR, sig_output_dir, n_cpus, 1, use_mp=True)
    signatures = find_signatures(sig_output_dir)

    print('Extracting potential regions.')
    regions = []
    for sig in signatures:
        i = np.where(info_df.index == sig)[0][0]
        start = max(0, i - offset)
        end = min(len(info_df), i + offset) + 1
        r = info_df.iloc[start:end]
        regions.append(r)
    
    return regions

def find_signatures(hmmsearch_output_dir):
    files = glob.glob(hmmsearch_output_dir + '/*.tab')
    signatures = set()

    for f in files:
        hmm_results = np.loadtxt(f, usecols=[0, 2, 4], dtype=np.str)

        if len(hmm_results) and len(hmm_results.shape) == 1:
            hmm_results = np.expand_dims(hmm_results, axis=0)
        
        signatures.update(prot for prot, _, _ in hmm_results)
    
    return sorted(signatures)

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-f', '--fasta', dest='fasta_file', help='Organism DNA Fasta file.', metavar='/path/to/organism.fa', type=str)
    parser.add_argument('-c', '--sequence-completeness', nargs='?', dest='sequence_completeness', help='Sequence completeness of DNA Fasta file. Available options: complete or partial.', metavar='seq_comp', choices=['complete', 'partial'], type=str)
    parser.add_argument('-o', '--output-directory', nargs='?', dest='output_dir', help='Output directory path.', metavar='output_dir', default='.', type=str)
    parser.add_argument('-n', '--number-of-cpus', nargs='?', dest='n_cpus', help='Number of CPUs to use.', default=1, type=int)
    parser.add_argument('-ho', '--hmmsearch-output-dir', nargs='?', dest='hmmsearch_output_dir', help='Hmmsearch output directory path.', metavar='hmmsearch_output_dir', default='./hmmsearch_output', type=str)
    args = parser.parse_args()

    if not os.path.exists(args.fasta_file):
        raise FileNotFoundError('No such file {}'.format(args.fasta_file))

    if not os.path.exists(CAS_HMM_DIR):
        extract_targz(CAS_HMM_TAR_PATH)
    
    if not os.path.exists(SIG_HMM_DIR):
        extract_targz(SIG_HMM_TAR_PATH)

    if not os.path.exists(GEN_HMM_DIR):
        extract_targz(GEN_HMM_TAR_PATH)

    if not os.path.exists(args.output_dir):
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    if not os.path.exists(args.hmmsearch_output_dir):
        Path(args.hmmsearch_output_dir).mkdir(parents=True, exist_ok=True)
    
    regions = find_potential_regions(args.fasta_file, args.sequence_completeness,
                                     args.output_dir, args.hmmsearch_output_dir, args.n_cpus)
