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

from pathlib import Path

from prodigal import run_prodigal, prodigal_fasta_to_genome_info
from utils import extract_targz

CAS_HMM_TAR_PATH = 'Cas_HMM.tar.gz'
GEN_HMM_TAR_PATH = 'Gen_HMM.tar.gz'

CAS_HMM_DIR = 'Cas_HMM'
GEN_HMM_DIR = 'Gen_HMM'

def preprocess_organism(fasta_file, sequence_completeness, output_dir):
    print('Running Prodigal on input Fasta file.')
    proteins_fasta_file = run_prodigal(fasta_file, sequence_completeness, output_dir)

    print('Generating Genome DataFrame')
    info_df = prodigal_fasta_to_genome_info(proteins_fasta_file, output_dir)

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-f', '--fasta', dest='fasta_file', help='Organism DNA Fasta file.', metavar='/path/to/organism.fa')
    parser.add_argument('-c', '--sequence-completeness', nargs='?', dest='sequence_completeness', help='Sequence completeness of DNA Fasta file. Available options: complete or partial.', metavar='seq_comp', choices=['complete', 'partial'])
    parser.add_argument('-o', '--output-directory', nargs='?', dest='output_dir', metavar='output_dir', default='.')
    args = parser.parse_args()

    if not os.path.exists(args.fasta_file):
        raise FileNotFoundError('No such file {}'.format(args.fasta_file))

    if not os.path.exists(CAS_HMM_DIR):
        extract_targz(CAS_HMM_TAR_PATH)

    if not os.path.exists(GEN_HMM_DIR):
        extract_targz(GEN_HMM_TAR_PATH)

    if not os.path.exists(args.output_dir):
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    preprocess_organism(args.fasta_file, args.sequence_completeness, args.output_dir)