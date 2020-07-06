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
import joblib
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
from wrapper import ClassifierWrapper

CAS_HMM_TAR_PATH = 'Cas_HMM.tar.gz'
SIG_HMM_TAR_PATH = 'Sig_HMM.tar.gz'
GEN_HMM_TAR_PATH = 'Gen_HMM.tar.gz'
ML_MODELS_TAR_PATH = 'ML_models.tar.gz'

CAS_HMM_DIR = 'Cas_HMM'
SIG_HMM_DIR = 'Sig_HMM'
GEN_HMM_DIR = 'Gen_HMM'
ML_MODELS_DIR = 'ML_models'

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
    output_dir += '/potential_regions/'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    regions = []
    for j, sig in enumerate(signatures):
        i = np.where(info_df.index == sig)[0][0]
        start = max(0, i - offset)
        end = min(len(info_df), i + offset) + 1
        r = info_df.iloc[start:end]
        r = r.assign(RefSignature=['no'] * r.shape[0])
        r.at[sig, 'RefSignature'] = 'yes'
        regions.append(r)
        r.to_csv(output_dir + f'/region_{j+1}.csv')
    
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

def find_boundaries(y_pred, max_gap, i_sig):
    last_left = i_sig
    left = i_sig - 1
    gap = 0

    while left >= 0 and gap <= max_gap:
        if y_pred[left] == 1:
            gap = 0
            last_left = left
        else:
            gap += 1
        
        left -= 1
    
    last_right = i_sig
    right = i_sig + 1
    gap = 0

    while right < len(y_pred) and gap <= max_gap:
        if y_pred[right] == 1:
            gap = 0
            last_right = right
        else:
            gap += 1
        
        right += 1
    
    return last_left, last_right

# def merge_overlapped_cassettes(cassettes_list, max_diff=2):
#     stop = False
#     merged_cassettes = []

#     while not stop:
#         stop = True

#         for i, cassette in enumerate(cassettes_list):
#             diff = np.array([len(set(cassette).symmetric_difference(set(other))) for other in cassettes_list], dtype=np.int)
#             indices_to_merge = np.where(diff <= max_diff)[0]

#             if len(indices_to_merge):
                            

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-f', '--fasta', dest='fasta_file', help='Organism DNA Fasta file.', metavar='/path/to/organism.fa', type=str)
    parser.add_argument('-c', '--sequence-completeness', nargs='?', dest='sequence_completeness', help='Sequence completeness of DNA Fasta file. Available options: complete or partial.', metavar='seq_comp', choices=['complete', 'partial'], type=str)
    parser.add_argument('-o', '--output-directory', nargs='?', dest='output_dir', help='Output directory path.', metavar='output_dir', default='.', type=str)
    parser.add_argument('-n', '--number-of-cpus', nargs='?', dest='n_cpus', help='Number of CPUs to use.', default=1, type=int)
    parser.add_argument('-g', '--maximum-gap', nargs='?', dest='max_gap', help='Maximum number of contiguous gaps allowed in a cassette (default: 1).', type=int, choices=range(6), default=1)
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
    
    if not os.path.exists(ML_MODELS_DIR):
        extract_targz(ML_MODELS_TAR_PATH)

    if not os.path.exists(args.output_dir):
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    if not os.path.exists(args.hmmsearch_output_dir):
        Path(args.hmmsearch_output_dir).mkdir(parents=True, exist_ok=True)
    
    proteins_fasta_file, regions_dataframes = find_potential_regions(args.fasta_file, args.sequence_completeness,
                                                          args.output_dir, args.hmmsearch_output_dir,
                                                          args.n_cpus)
    
    regions_fasta_files = extract_regions_sequences(proteins_fasta_file, regions_dataframes, args.output_dir)
    gen_output_dir = args.hmmsearch_output_dir + '/' + GEN_HMM_DIR
    gen_hmm_names = np.loadtxt(GEN_HMM_NAMES_FILE, dtype=np.str)
    gen_hmm_to_index = dict(zip(gen_hmm_names, np.arange(len(gen_hmm_names))))
    regions_hmm_list = []
    regions_prop_list = []
    
    print('Extracting Generic HMM features and protein properties features.')
    for r, f in zip(regions_dataframes, regions_fasta_files):
        r_name = f.rsplit('/', 1)[1].replace('.fasta', '')
        run_hmmsearch(f, GEN_HMM_DIR, gen_output_dir + '/' + r_name, args.n_cpus, 1000, use_mp=False)

        X_hmm = build_hmm_matrix(r, r_name, args.output_dir + '/potential_regions',
                                 args.hmmsearch_output_dir + '/Gen_HMM/' + r_name, gen_hmm_to_index)
        
        regions_hmm_list.append(X_hmm)

        X_prop = build_protein_properties_matrix(r, f, r_name, args.output_dir + '/potential_regions')

        regions_prop_list.append(X_prop)
    
    print('Finding cassette boundaries.')
    scaler_ert = joblib.load(ML_MODELS_DIR + '/ert_bound/MaxAbsScaler_hmm.joblib')
    svd_ert = joblib.load(ML_MODELS_DIR + '/ert_bound/TruncatedSVD_hmm.joblib')
    clf_ert = joblib.load(ML_MODELS_DIR + '/ert_bound/ExtraTreesClassifier_hmm.joblib')

    scaler_dnn = joblib.load(ML_MODELS_DIR + '/dnn_bound/MaxAbsScaler_hmm.joblib')
    svd_dnn = joblib.load(ML_MODELS_DIR + '/dnn_bound/TruncatedSVD_hmm.joblib')
    clf_dnn = joblib.load(ML_MODELS_DIR + '/dnn_bound/MLPClassifier_hmm.joblib')

    ert_boundaries_list = []
    dnn_boundaries_list = []

    for rdf, X_hmm in zip(regions_dataframes, regions_hmm_list):
        i_signature = np.where(rdf['RefSignature'] == 'yes')[0][0]
        i_signature_rep = [i_signature] * X_hmm.shape[0]
        X_sig = X_hmm[i_signature_rep]
        X_hmm = sparse.hstack((X_sig, X_hmm))

        X_hmm_ert = scaler_ert.transform(X_hmm)
        X_hmm_ert = svd_ert.transform(X_hmm_ert)
        y_pred_ert = clf_ert.predict(X_hmm_ert)
        ert_boundaries = find_boundaries(y_pred_ert, args.max_gap, i_signature)
        ert_boundaries_list.append(ert_boundaries)
        print(ert_boundaries)

        X_hmm_dnn = scaler_dnn.transform(X_hmm)
        X_hmm_dnn = svd_dnn.transform(X_hmm_dnn)
        y_pred_dnn = clf_dnn.predict(X_hmm_dnn)
        dnn_boundaries = find_boundaries(y_pred_dnn, args.max_gap, i_signature)
        dnn_boundaries_list.append(dnn_boundaries)
        print(dnn_boundaries)
        print()

    # print('Labeling Cas proteins.')
    # maxabs_scaler_ert = joblib.load(ML_MODELS_DIR + '/ert_prot/MaxAbsScaler_hmm2019_prop.joblib')
    # minmax_scaler_ert = joblib.load(ML_MODELS_DIR + '/ert_prot/MinMaxScaler_hmm2019_prop.joblib')
    # clf_ert = joblib.load(ML_MODELS_DIR + '/ert_prot/ClassifierWrapper_hmm2019_prop.joblib')

    # maxabs_scaler_dnn = joblib.load(ML_MODELS_DIR + '/dnn_prot/MaxAbsScaler_hmm2019_prop.joblib')
    # minmax_scaler_dnn = joblib.load(ML_MODELS_DIR + '/dnn_prot/MinMaxScaler_hmm2019_prop.joblib')
    # clf_dnn = joblib.load(ML_MODELS_DIR + '/dnn_prot/ClassifierWrapper_hmm2019_prop.joblib')

    # final_predictions_ert = []
    # final_predictions_dnn = []

    # for i in range(len(regions_dataframes)):
    #     X_hmm = regions_hmm_list[i].toarray()
    #     X_prop = regions_prop_list[i]

    #     start_ert, end_ert = ert_boundaries_list[i]
    #     X_hmm_ert = X_hmm[start_ert:end_ert+1]
    #     X_hmm_ert = maxabs_scaler_ert.transform(X_hmm_ert)
    #     X_prop_ert = X_prop[start_ert:end_ert+1]
    #     X_prop_ert = minmax_scaler_ert.transform(X_prop_ert)
    #     X_ert = np.hstack((X_hmm_ert, X_prop_ert))
    #     y_pred_ert = clf_ert.predict(X_ert)

    #     start_dnn, end_dnn = dnn_boundaries_list[i]
    #     X_hmm_dnn = X_hmm[start_dnn:end_dnn+1]
    #     X_hmm_dnn = maxabs_scaler_dnn.transform(X_hmm_dnn)
    #     X_prop_dnn = X_prop[start_dnn:end_dnn+1]
    #     X_prop_dnn = minmax_scaler_dnn.transform(X_prop)
    #     X_dnn = np.hstack((X_hmm_dnn, X_prop_dnn))
    #     y_pred_dnn = clf_dnn.predict(X_dnn)
    
