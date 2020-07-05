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
from Bio import SeqIO
from scipy import sparse

# Project imports
from prodigal import run_prodigal, prodigal_fasta_to_genome_info
from hmmer import run_hmmsearch
from utils import extract_targz
from protein_features import get_protein_features

CAS_HMM_TAR_PATH = 'Cas_HMM.tar.gz'
SIG_HMM_TAR_PATH = 'Sig_HMM.tar.gz'
GEN_HMM_TAR_PATH = 'Gen_HMM.tar.gz'

CAS_HMM_DIR = 'Cas_HMM'
SIG_HMM_DIR = 'Sig_HMM'
GEN_HMM_DIR = 'Gen_HMM'

GEN_HMM_NAMES_FILE = 'gen_hmm_names.txt'

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
        r = r.assign(RefSignature=['no'] * r.shape[0])
        r.at[sig, 'RefSignature'] = 'yes'
        regions.append(r)
    
    return proteins_fasta_file, regions

def find_signatures(hmmsearch_output_dir):
    files = glob.glob(hmmsearch_output_dir + '/*.tab')
    signatures = set()

    for f in files:
        hmm_results = np.loadtxt(f, usecols=[0, 2, 4], dtype=np.str)

        if len(hmm_results) and len(hmm_results.shape) == 1:
            hmm_results = np.expand_dims(hmm_results, axis=0)
        
        signatures.update(prot for prot, _, _ in hmm_results)
    
    return sorted(signatures)

def extract_regions_sequences(proteins_fasta_file, regions, output_dir):
    output_dir = output_dir + '/potential_regions'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    regions_fasta_files = []
    
    with open(proteins_fasta_file, 'r') as file_:
        sequences = SeqIO.to_dict(SeqIO.parse(file_, 'fasta'))
        
    for i, r in enumerate(regions):
        r_seqs = [sequences[prot_id] for prot_id in r.index]

        out_path = output_dir + f'/region_{i+1}.fasta'
        with open(out_path, 'w') as out_file:
            for s in r_seqs:
                out_file.write(s.format('fasta'))                         
        
        regions_fasta_files.append(out_path)
    
    return regions_fasta_files

def build_hmm_matrix(region, region_name, matrix_output_dir, hmmsearch_output_dir, hmm_to_index):
    prot_to_index = dict(zip(region.index, np.arange(region.shape[0])))
    hmm_results_files = glob.glob(hmmsearch_output_dir + '/*.tab')
    X_hmm = sparse.dok_matrix((region.shape[0], len(hmm_to_index)), dtype=np.double)
            
    for hmm_f in hmm_results_files:
        hmm = np.loadtxt(hmm_f, usecols=[0, 2, 5], dtype=np.str)
        
        if len(hmm) and len(hmm.shape) == 1:
            hmm = np.expand_dims(hmm, axis=0)

        for prot, name, bitscore in hmm:
            i = prot_to_index[prot]
            j = hmm_to_index[name]
            bitscore = float(bitscore)
            X_hmm[i, j] = max(X_hmm[i, j], bitscore)
        
    X_hmm = sparse.csr_matrix(X_hmm)
    sparse.save_npz(matrix_output_dir + '/' + region_name + '.npz', X_hmm, compressed=True)
    
    return X_hmm

def build_protein_properties_matrix(region, region_fasta, region_name, matrix_output_dir):
    with open(region_fasta, 'r') as file_:
        sequences = SeqIO.to_dict(SeqIO.parse(file_, 'fasta'))
    
    X_prop = []
    for prot_id in region.index:
        seq_record = sequences[prot_id]
        names, prot_features = get_protein_features(str(seq_record.seq))
        X_prop.append(prot_features)
    
    X_prop = np.array(X_prop)
    np.savetxt(matrix_output_dir + '/' + region_name + '.prop.txt', X_prop, header=' '.join(names))
    
    return X_prop

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
    
    proteins_fasta_file, regions = find_potential_regions(args.fasta_file, args.sequence_completeness,
                                                          args.output_dir, args.hmmsearch_output_dir,
                                                          args.n_cpus)
    
    regions_fasta_files = extract_regions_sequences(proteins_fasta_file, regions, args.output_dir)
    gen_output_dir = args.hmmsearch_output_dir + '/' + GEN_HMM_DIR
    gen_hmm_names = np.loadtxt(GEN_HMM_NAMES_FILE, dtype=np.str)
    gen_hmm_to_index = dict(zip(gen_hmm_names, np.arange(len(gen_hmm_names))))
    regions_hmm_list = []
    regions_prop_list = []
    
    print('Extracting Generic HMM features and protein properties features.')
    for r, f in zip(regions, regions_fasta_files):
        r_name = f.rsplit('/', 1)[1].replace('.fasta', '')
        run_hmmsearch(f, GEN_HMM_DIR, region_output_dir, args.n_cpus, 1000, use_mp=False)

        X_hmm = build_hmm_matrix(r, r_name, args.output_dir + '/potential_regions',
                                 args.hmmsearch_output_dir, gen_hmm_to_index)
        
        regions_hmm_list.append(X_hmm)

        X_prop = build_protein_properties_matrix(r, f, r_name, args.output_dir + '/potential_regions')

        regions_prop_list.append(X_prop)
    
      

