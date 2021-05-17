from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import PDB
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
import pickle
import os
from os import path
import sys
import warnings
import numpy as np
import numpy
import subprocess
import collections
import pandas as pd
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
if not sys.warnoptions:
    warnings.simplefilter("ignore")

from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from imblearn.combine import SMOTETomek
from imblearn.over_sampling import RandomOverSampler
from sklearn.model_selection import cross_validate
from copy import deepcopy

def calc_residue_dist(residue_one, residue_two) :
    """Computes and returns the distance between two residues, by comparing the position of their alpha carbons"""
    #TODO : return an integer representing the distance between the two residues, in Angstrom
    try:
        distance  = residue_one["CA"] - residue_two["CA"]
    except:
        distance = 0
    return distance

def compute_distance_matrix(residues) :
    """Computes a matrix of size len(Seq) * len(Seq) with distances between each residue."""

    #TODO : return a numpy 2D array of distances between each residue of the structure.
    #Tip : you might want to make sure you filter out consecutive residues at this step.
    distance_matrix = []
    for residue1 in residues:
        distance_line = []
        for residue2 in residues:
            distance = abs(residue1["CA"] - residue2["CA"])
            distance_line.append(distance)
        distance_matrix.append(distance_line)
    return distance_matrix

def extract_residues(model):
    """Returns a list of protein residues given a PDB model"""

    #TODO : return a list of protein residues given a PDB model
    residues = []
    for chain in model:
        for residue in chain:
            if PDB.is_aa(residue, standard=True): residues.append(residue)
    return residues

def get_dssp_info(PDB_file,model,dir):
    """Runs DSSP on protein input"""
    #TODO : you can run DSSP through biopython. The output contains a lot of useful information.
    #Tip : make sure your secondary structure indexing matches the sequence order in the PDB file!
    dssp = dssp_dict_from_pdb_file(dir+'/'+PDB_file)
    return dssp[0]

def write_fasta(sequence,PDB_file):
    """Writes sequence to fasta file named after the PDB the sequence comes from"""

   #TODO : implement the writing of a fasta file from the sequence obtained from the PDB file.

   #NOTE: writes all sequences to a fasta file called "all.fasta" -> use this combined fasta to run tmhmm
    filename = os.getcwd() + '/data/all.fasta'
    f = open(filename, "a")
    f.write(">"+PDB_file[:-4]+'\n'+str(sequence)+'\n')
    f.close()
    return PDB_file[:-4]

# def run_tmhmm(pdb_file, filename):
#     """Runs tmhmm on input fasta files, sends the output to STDout, processes output into secondary structure"""
#     #TODO : run TMHMM and output its resulting secondary structure.
#     #If you cannot get tmhmm to run, you can write all your sequences to a FASTA file,
#     # use the webserver and parse the output file with this function. Otherwise, you can use this function to run
#     #TMHMM on some FASTA file and parse its output.

#     #NOTE: ran tmhmm on all.fasta using webserver
#     # in this function, parsing through output file (called tmhmm_results.txt)
#     df = pd.read_csv(filename, sep='\t', header=None)
#     col = df.loc[df[0] == pdb_file]
#     length = int(col[1].reset_index(drop=True)[0][4:])
#     idx = col[5].reset_index(drop=True)[0][9:]
#     idx = re.split('o|i|-',idx)
#     # NOTE: first position = 1 in TMHMM (NOT 0)
#     idx.pop(0)
#     idx.pop(-1)
#     ss = ''
#     prev = 0
#     num = 0
#     if idx:
#         for i in idx:
#             num+=1
#             if num%2 == 0:
#                 ss+=('H'*(int(i)-prev+1))
#             else:
#                 ss+=('C'*(int(i)-prev-1))
#             prev = int(i)
#         # last index until end
#         ss+=('C'*(length-int(idx[-1])))
#     else:
#         ss+=('C'*length)
#     return ss

def get_contact_numbers(contact_map):
    """Returns the proportion of residues involved in intramolecular contacts"""
    #TODO : a utility function for the number of contacts in the protein, which you might use to make sure your output makes sense
    return numpy.count_nonzero(contact_map) / float(contact_map.size)

#HERE
def generate_ML_dataset(sequence, hydro, dssp_ss, dssp_asa, dssp_phi, dssp_psi,has_contact,DSSP_vector, oracle):
    # def generate_ML_dataset(sequence, hydro, dssp_ss, dssp_asa, dssp_phi, dssp_psi, tm_ss,has_contact,DSSP_vector, TMHMM_vector, oracle):
    """generates vectors for machine learning based on sequence, DSSP and TMHMM outputs"""

    ########################################################################################
    #TODO : generate machine learning dataset from PDB info, and append it to the three vectors in input.

    #This function takes as input the sequence, the two secondary structures, and the has_contact boolean array.
    #use the first 4 features to append elements defined as follows to the two vectors and the oracle.

    # the features are a list of 9 tuples ("AA","SS of AA"). The oracle states whether there is intramolecular contact (<5 Ang) at the center of this subsequence of length 9.
    #An example vector element is the following : [("Y","H),("S","C"),("A","C"),("S","C"),("A","C"),("S","C"),("A","C"),("Y","H)]
    #The ith list of tuples of DSSP_vector and TMHMM_vector should describe the same sequence.

    # the oracle is a list of booleans values, establishing presence of a contact for each vector.
    # the intramolecular contact label of the ith list of tuples in DSSP_vector and TMHMM is found at index i of the oracle.

    #return the three vectors after appending the new elements.
    ########################################################################################
    sequence = str(sequence)
    seq_nine = list(split_by_n(sequence, 9))
    # tm_nine = list(split_by_n(tm_ss, 9))
    dssp_nine = list(split_by_n(dssp_ss, 9))
    dssp_asa_nine = list(split_by_n(dssp_asa,9))
    dssp_phi_nine = list(split_by_n(dssp_phi,9))
    dssp_psi_nine = list(split_by_n(dssp_psi,9))
    hydro_nine = list(split_by_n(hydro,9))

    for i in enumerate(seq_nine):
        if(len(seq_nine[i[0]]) < 9): continue
        DSSP_vector.append(zip(seq_nine[i[0]], dssp_nine[i[0]], dssp_asa_nine[i[0]], dssp_phi_nine[i[0]], dssp_psi_nine[i[0]], hydro_nine[i[0]]))
        # TMHMM_vector.append(zip(seq_nine[i[0]], tm_nine[i[0]]))

    i = 0
    while i<len(sequence):
        if i == 0:
            oracle.append(has_contact[i+4])
        else:
            oracle.append(has_contact[i])
        i = i + 9
    # return DSSP_vector, TMHMM_vector, oracle
    return DSSP_vector, oracle

def split_by_n(seq, n):
    '''Divides a sequence into n subunits'''
    while seq:
        yield seq[:n]
        seq = seq[n:]

def removeConsecutives(matrix):
    contact_map = []
    row_idx = 0
    for row in matrix:
        contact_list = []
        col_idx = 0
        for column in row:
            if (column < 5 and abs(col_idx - row_idx) >= 10):
                contact_list.append(True)
            else:
                contact_list.append(False)
            col_idx += 1
        row_idx += 1
        contact_map.append(contact_list)
    contact_map = np.array(contact_map)
    return(contact_map)

#HEFE
def get_PDB_info(dir):
    """Extracts sequence, DSSP secondary structure, TMHMM secondary structure and contact information from PDB files in input directory"""

    #the three vectors you are required to fill.
    # DSSP_vector, TMHMM_vector, oracle = [],[],[]
    DSSP_vector, oracle = [],[]

    print("There are",len(os.listdir(dir)),"PDB files to parse")
    if os.path.exists(os.getcwd()+'/data/all.fasta'):
        os.remove(os.getcwd()+'/data/all.fasta')

    #Assemble a machine learning dataset incrementally, for each PDB file in the directory
    for ind,PDB_file in enumerate(os.listdir(dir)):
        if ind%10==0:
            print("Working on structure",ind)
        print(PDB_file)
        #Step 1 : parse your PDB file with biopython to obtain a model object
        p = PDB.FastMMCIFParser()
        path = os.path.join(dir, PDB_file)
        structure = p.get_structure(PDB_file[:-4].upper(), path)
        model = structure[0]


        #TODO : extract a list of residues from your model object
        residues = extract_residues(model)
        matrix = compute_distance_matrix(residues)
        contact_map = removeConsecutives(matrix)
        has_contact = [True if True in contact_map[residue] else False for residue in contact_map]

        #TODO : contact info should return the proportion of residues that have an intramolecular contact in your object.
        contact_info = get_contact_numbers(contact_map)
        # print(contact_info,"contacts")

        #TODO : obtain the secondary structure prediction of the PDB model with DSSP
        dssp_info = get_dssp_info(PDB_file,model,dir)

        #TODO : obtain the sequence of the PDB file in some way of your choice.
        sequence = ""
        dssp_ss = "" #ss stands for secondary structure

        # get sequence from structure (called model)
        ppb = PDB.PPBuilder()
        for pp in ppb.build_peptides(model):
            sequence+=(pp.get_sequence())

        # added hydrophobicity (uses Kyte-Doolittle scale)
        kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
           'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
           'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
           'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
        hydro = []
        for residue in sequence:
            hydro.append(kd[residue])
        # get dssp secondary structure (if using biopython)
        #NEW: from dssp, also get relative ASA
        dssp_asa = []
        dssp_phi = []
        dssp_psi = []
        for aa in sorted(dssp_info.keys()):
            ss = dssp_info[aa][1]
            #add relative ASA
            asa = dssp_info[aa][2]
            dssp_asa.append(asa)
            #add phi and psi
            dssp_phi.append(dssp_info[aa][3])
            dssp_psi.append(dssp_info[aa][4])
            # H = G is 3-turn, H = 4-turn, I = 5-turn helix
            if (ss in ['H', 'G', 'I']):
                dssp_ss+="H"
            #everything that isnt a helix is considered a coil...
            elif (ss == 'E'):
                dssp_ss+="E"
            else:
                dssp_ss+="C"

        #TODO : write the sequence to a fasta file to call TMHMM with it, or to use the webserver
        # filename = write_fasta(sequence,PDB_file)

        # already ran tmhmm on compiled fasta file on webserver -> need to parse file
        #TODO : obtain secondary structure prediction for this FASTA file with TMHMM
        # tm_ss = run_tmhmm(filename, os.getcwd()+'/tmhmm_results.txt')

        # DSSP_vector, TMHMM_vector, oracle = generate_ML_dataset(sequence, hydro, dssp_ss, dssp_asa, dssp_phi, dssp_psi, tm_ss,has_contact,DSSP_vector, TMHMM_vector, oracle)
        DSSP_vector, oracle = generate_ML_dataset(sequence, hydro, dssp_ss, dssp_asa, dssp_phi, dssp_psi,has_contact,DSSP_vector, oracle)

    # return DSSP_vector, TMHMM_vector, oracle
    return DSSP_vector, oracle

def generate_dataset():
    """Runs the PDB parsing utility"""
    # DSSP_vector, TMHMM_vector, oracle = get_PDB_info(os.path.join("data","pdb"))
    DSSP_vector, oracle = get_PDB_info(os.path.join("data","pdb"))
    print("successfully parsed PDB data")

    #store a pickle of your results to avoid repeating get_PDB_info
    # pickle.dump((DSSP_vector, TMHMM_vector, oracle),open("ML_ready_dataset.pickle","wb"))
    # return DSSP_vector, TMHMM_vector, oracle
    pickle.dump((DSSP_vector, oracle),open("ML_ready_dataset.pickle","wb"))
    return DSSP_vector, oracle

def split_dataset(X, Y):
    """Splits the dataset into training and testing"""
    # TODO : split the X and Y dataset into a reasonable training set and a test set. Your test set should have have 20% of the datapoints.
    # Tip : look up train_test_split with scikit-learn
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)
    return X_train, X_test, Y_train, Y_test

# def format_simple_dataset(vector,solutions):
def format_simple_dataset(alphads, betads, aandbds, adashbds, gpcrds):
    """takes as input a vector of sequence strings and a vector of booleans
    outputs a vector of size 9 vectors of tuples and an array of binary numbers"""

    lst = [alphads, betads, aandbds, adashbds, gpcrds]
    xlst, ylst = [],[]
    for ds in lst:
        vector_orig, solutions = ds
        vector = deepcopy(vector_orig)
        AA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T", "X"]
        contact_truth = [False, True]

        training_size = len(vector)
        formatted_X = []
        for i in vector[:training_size]:
            current_vec = np.zeros(189)
            for ind, j in enumerate(i):
                nuc = AA.index(j[0])
                current_vec[21 * ind + nuc] = 1
            current_vec = np.array(current_vec)
            formatted_X.append(current_vec)
        xlst.append(np.array(formatted_X))
        formatted_Y = np.zeros(training_size)
        for i in range(training_size):
            formatted_Y[i] = contact_truth.index(solutions[i])
        ylst.append(np.array(formatted_Y))
    return xlst, ylst

# def format_one_hot_dataset(vector,solutions):
def format_one_hot_dataset(alphads, betads, aandbds, adashbds, gpcrds):
    """Takes as input a vector of 9 sequence, ss tuples, outputs a vector of one-hot vectors of size 198 and a binary vector"""
    AA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T", "X"]
    ss = ["H", "C", "E"]
    contact_truth = [False, True]

    lst = [alphads, betads, aandbds, adashbds, gpcrds]
    xlst, ylst = [],[]
    for ds in lst:
        vector, solutions = ds
        training_size = len(vector)
        formatted_X = []
        for i in vector[:training_size]:
            current_vec = np.zeros(198)
            for ind, j in enumerate(i):
                nuc = AA.index(j[0])
                sec = ss.index(j[1])
                current_vec[22 * ind + nuc] = 1
                current_vec[22 * ind + 21] = sec
            current_vec = np.array(current_vec)
            formatted_X.append(current_vec)
        xlst.append(np.array(formatted_X))
        formatted_Y = np.zeros(training_size)
        for i in range(training_size):
            formatted_Y[i] = contact_truth.index(solutions[i])
        ylst.append(np.array(formatted_Y))
    return xlst, ylst

# def run_NN_on_sequence(vector, solutions):
def run_NN_on_sequence(alphads, betads, aandbds, adashbds, gpcrds):
    """Trains a scikit-learn basic 4 layers neural network on a sequence vector, outputs results"""
    #X is a the input vector of ML-ready numpy arrays, Y is the corresponding oracle
    Xlst, Ylst = format_simple_dataset(alphads, betads, aandbds, adashbds, gpcrds)
    #TODO : fill in split_dataset


    X_train_total, Y_train_total = [], []
    X_test_lst, Y_test_lst = [], []
    for i in range(0, len(Xlst)):

        #currently using SMOTETomek
        method = SMOTETomek()
        X_resampled, y_resampled = method.fit_resample(Xlst[i], Ylst[i])
        X_train, X_test, Y_train, Y_test = split_dataset(X_resampled, y_resampled)

        # oversample = RandomOverSampler(sampling_strategy='minority')
        # X_resampled, y_resampled = oversample.fit_resample(Xlst[i], Ylst[i])
        # X_train, X_test, Y_train, Y_test = split_dataset(X_resampled, y_resampled)

        X_train_total.append(X_train)
        Y_train_total.extend(Y_train)
        X_test_lst.append(X_train)
        Y_test_lst.append(Y_train)
    X_train = np.vstack(X_train_total)
    Y_train = np.vstack(Y_train_total)

    clf = MLPClassifier(solver='sgd', alpha=1e-6, hidden_layer_sizes=(156), activation="logistic", shuffle=True,
                        verbose=False, random_state=1, tol=1e-5, max_iter=350)
    clf.fit(X_train, Y_train)

    print('done fit')
    print("score, sequence only")
    for i in range(0, 5):
        Y_pred_train = clf.predict(X_train)
        Y_pred_test = clf.predict(X_test_lst[i])
        if i == 0:
            print("ALPHA ONLY")
        elif i == 1:
            print("BETA ONLY")
        elif i == 2:
            print("A+B")
        elif i == 3:
            print("A/B ")
        else:
            print("GPCR")
        print("TEST SCORES", classification_report(Y_test_lst[i], Y_pred_test))
        print("CONFUSION MATRIX\n", confusion_matrix(Y_test_lst[i], Y_pred_test))

# def run_NN_for_ss_type(vector,solutions):
def run_NN_for_ss_type(alphads, betads, aandbds, adashbds, gpcrds):
    """Trains a scikit-learn basic 4 layers neural network on a sequence-structure vector, outputs results"""
    #X is a the input vector of ML-ready numpy arrays, Y is the corresponding oracle
    Xlst, Ylst = format_one_hot_dataset(alphads, betads, aandbds, adashbds, gpcrds)
    #TODO : fill in split_dataset

    X_train_total, Y_train_total = [], []
    X_test_lst, Y_test_lst = [], []
    for i in range(0, len(Xlst)):
        #using SmoteTOMEK
        method = SMOTETomek()
        X_resampled, y_resampled = method.fit_resample(Xlst[i], Ylst[i])
        X_train, X_test, Y_train, Y_test = split_dataset(X_resampled, y_resampled)

        # oversample = RandomOverSampler(sampling_strategy='minority')
        # X_resampled, y_resampled = oversample.fit_resample(Xlst[i], Ylst[i])
        # X_train, X_test, Y_train, Y_test = split_dataset(X_resampled, y_resampled)
        
        X_train_total.append(X_train)
        Y_train_total.extend(Y_train)
        X_test_lst.append(X_train)
        Y_test_lst.append(Y_train)
    X_train = np.vstack(X_train_total)
    Y_train = np.vstack(Y_train_total)

    clf = MLPClassifier(solver='sgd', alpha=1e-6, hidden_layer_sizes=(156), activation="logistic", shuffle=True,
                        verbose=False, random_state=1, tol=1e-5, max_iter=350)
    clf.fit(X_train, Y_train)

    print('done fit')
    print("score, all features")
    for i in range(0, 5):
        Y_pred_train = clf.predict(X_train)
        Y_pred_test = clf.predict(X_test_lst[i])
        if i == 0:
            print("ALPHA ONLY")
        elif i == 1:
            print("BETA ONLY")
        elif i == 2:
            print("A+B")
        elif i == 3:
            print("A/B ")
        else:
            print("GPCR")
        print("TEST SCORES", classification_report(Y_test_lst[i], Y_pred_test))
        print("CONFUSION MATRIX\n", confusion_matrix(Y_test_lst[i], Y_pred_test))

def predict_intramolecular_contacts(alphads, betads, aandbds, adashbds, gpcrds):
    """Compares neural network results for DSSP and TMHMM secondary structures"""
    print("SEQUENCE ONLY RESULTS")
    #use either vector for sequence only prediction
    run_NN_on_sequence(alphads, betads, aandbds, adashbds, gpcrds)
    print("DSSP RESULTS")
    run_NN_for_ss_type(alphads, betads, aandbds, adashbds, gpcrds)

def get_data():
    #GPCR
    L = ["1D6G", "1FJR", "1MF6", "1NZS", "1WSO", "1Y7J", "1YM7", "1YTV", "2DCO", "2JX0", "2JX4", "2JX9", "2K3U", "2K9P", "2KI9", "2KOE", "2L60", "2L63", "2LMK", "2LNL", "2LQ4", "2M9O", "2N2F", "2N55", "2NZ1", "2O8Z", "2QKH", "2RGN", "2RH1", "2RRS", "2X57", "2XVT", "2XWT", "2Y00", "2YDV", "2YX8", "3EML", "3H3G", "3I8S", "3N7S", "3N93", "3NKQ", "3ODU", "3PBL", "3RZE",
"3T33", "3UGU", "3UON", "3V2Y", "4DJH", "4DLO", "4DLQ", "4EIY", "4I6O", "4IAR", "4IB4", "4JKV", "4JQI", "4K5Y", "4L6R", "4MBS", "4MQW", "4N6H", "4NUU", "4OR2", "4P39", "4P3A", "4PNI", "4PXZ", "4QIM", "4R7V", "4R7X", "4RMK", "4TND", "4U14", "4U15", "4X1H", "4XES", "4XNV", "4XT1", "4Z35", "4ZUD", "5CXV", "5DHG", "5DSG", "5EE7", "5FTT", "5KVM", "5KZZ", "5NDD", "5NX2", "5O9H", "5OLL", "5OTU", "5T04", "5T1A", "5UEN", "5UNF", "5V57", "5VBL", "5VEW",
"5WB1", "5WIU", "5XJM", "5ZKC", "5ZTY", "6AK3", "6B73", "6BD4", "6BQH", "6C1R", "6CM4", "6CRK", "6D1U", "6D27", "6D32", "6DO1", "6F27", "6F3V", "6F3W", "6F3X", "6F3Y", "6FFI", "6FJ3", "6GPS", "6GPX", "6HLP", "6IDX", "6IIU", "6J20", "6M9T"]
    
    #additional data
    # L=['6MHQ', '6MHY', '2ZW3', '5VOT', '6L3U', '6DM0', '5GJW', '6AKF', '5B2G', '6QKZ', '3X29', '6OV2', '4P79', '2BL8', '1TDP', '3KAW', '5FIG', '5FJE', '6EI0', '3D19', '3DBY', '6EK9', '3LMF', '2C12', '1EGD', '3MDE', '2D29', '1IVH', '1SIQ', '1RX0', '1U8V', '1IS2', '1BUC', '1JQI', '1UKW', '1W07', '1R2J', '2CJ4', '1X91', '1XG2', '1PI1', '5B6B', '5NCM', '5NCN', '1R3B', '2ETD', '1V9V', '2M9X', '2FEF', '2FUG', '2QZG', '2QSB', '2HI7', '2GSC', '2RLD', '2HGK', '1GKZ', '1Y8O', '1JM6', '3RCW', '5UVW', '3G0L', '3DWY', '3MB3', '4UYG', '6ULQ', '4TU6', '4YYG', '5FFV', '2YYN', '1WUM', '3I3J', '5H6Y', '1EQF', '5A7C', '4QBM', '2I7K', '4QY4', '5EA1', '2RFJ', '3D7C', '5K29', '5V84', '2DKW', '1EQF', '4QNS', '5UIY', '5VS7', '6NIM', '4PY6', '6LQX', '2DWW', '5W0I', '5N13', '6DF7', '4PKL', '4LDF', '2WP2', '5IGL', '6CW0', '3FKM', '3TLP', '5C4Q', '6MF9', '3MQM', '5C8G', '5K04', '5FEA', '6SW3', '5TCK', '1E6I', '5ML0', '2L5E', '3Q23', '2YQD', '5CQZ', '6NEY', '4L3U', '5CD4', '3WVO', '4H3T', '2WL8', '2HJM']
    # L = ['2B5L', '4D5L', '4PC0', '1TBG', '1K8K', '5IFE', '2YNP', '2CE8', '1YFQ', '6FYX', '1VYH', '1SQ9', '1P22', '1NEX', '5LQW', '2YBA', '6TBL', '2OVR', '3DWL', '2XL2', '4J87', '4BTS', '1PGU', '1ERJ', '3CFV', '4PSW', '3MKQ', '6K42', '1NR0', '5HY7', '4XYH', '2AQ5', '3EU7', '6UHC', '3J7R', '4V3P', '4V5O', '5OPT', '4OZU', '2YNP', '6YW7', '4V6W', '7KYX', '2BBK', '2MAD', '1L0Q', '4R5O', '1JMX', '1PBY', '5A9Q', '4XMM', '2OIT', '4FHN', '5A9Q', '1XIP', '4Q9T', '1XKS', '5HAX', '4MHC', '2XBG', '2EBS', '5MUK', '5OJ3', '4IRT', '5MUM', '5MVH', '4D9S', '1A12', '1JTD', '3MVD', '2HU7', '1QFM', '2BWR', '4O9X', '4TQM', '4IGL', '5MB4', '6NAU', '3FGB', '4QRJ', '1RI6', '1SHY', '3OKY', '1Q47', '6FKK', '1OLZ', '6QP9', '6WTS', '1TYE', '1M1X', '4WJK', '1K32', '1JOF', '1GOF', '1QNI', '1FWX', '1UTC', '5MXF', '6RG2', '7BB4', '5C9O', '6T96', '3S25', '4OU9', '4RSC', '5U8Y', '6VCF', '6N20', '6OJT', '5V2D', '5J53', '3NPE', '4I0O']
    # L = ['6LFB', '1ES9', '1ZMB', '1K7C', '3BZW', '1FLC', '4NRD', '1Z8H', '2HSJ', '2APJ', '7C23', '1ESC', '4M8K', '4JGG', '4O8V', '4KNC', '7C2A', '7C2C', '7C2D', '5A4A', '2O14', '1YZF', '1FXW', '3DC7', '4I8I', '3HP4', '2QWX', '1D4A', '1JA1', '1BU5', '1NNI', '3S2Y', '2RG1', '5NLL', '1YDG', '1CZN', '1YCG', '3W78', '4N9Q', '1B1C', '2Z9D', '1VME', '2Q62', '2ZKI', '5MP4', '6JXN', '2ARK', '4M0C', '3N3A', '4LAF', '2A5L', '1RLI', '2FZV', '1SQS', '1AG9', '1DXQ', '3HR4', '6DXP', '2HPV', '1YOB', '1RTT', '5F51', '2BMX', '5JRO', '4C0X', '4BMP', '1T0I', '1T5B', '1TLL', '1RCF', '1E5D', '1QRD', '3SVL', '4DY4', '2MT9', '4GI5', '2M6R', '3W77', '3U7R', '6QU0', '6OHK', '2XOE', '6EBQ', '1YKG', '1FUE', '1RLJ', '2FZ5', '2FCR', '1BVY', '2WC1', '3POR', '4BMO', '3D7N', '7BX9']
    # L = ['3ZNJ', '2CB2', '5XZQ', '3LO3', '1TR0', '6FKS', '4AU9', '3VEF', '5LOQ', '2CYY', '5B09', '2ZDP', '3HX9', '3HDS', '1MLI', '1LQ9', '6VSA', '4HL9', '1VQY', '2OMO', '2D3Q', '2AP6', '2CFX', '5KAK', '3KG0', '2D3Q', '2HIQ', '4LBH', '1T0T', '1VDH', '5T2K', '1T0T', '1VDH', '1X8D', '2PD1', '1Q4R', '1XBW', '4ZOS', '2QLW', '2ZDO', '3QMQ', '4GY7', '5DE0', '5GT2', '5VJ0', '6FSK', '5C2I', '3E8O', '4GU7', '5DE0', '5GT2', 'FVJ0', '6FSK', '5C2I', '1TZ0', '1X7V', '2CG4', '2IIZ', '2IFX', '2OKQ', '1IUJ', '1MWQ', '1RJJ', '2FIU', '2FTR', '1TUV', '2IIZ', '1Y0H', '1I1G', '2JSX', '4DPO', '2GFF', '4G2C', '3BGU', '2QYC', '4DN9', '4G2C', '1SQE', '7BIO', '2GVK', '2GO8', '1S7I', '1Q8B', '1TUW2GVK', '5K9F', '5IXU', '2BBE', '2FB0', '3BXV', '3BM7', '6HHN']
    pdbl = PDB.PDBList(verbose=False)

    for pdb in L:
        path_name = pdbl.retrieve_pdb_file(pdb, pdir='data/pdb/', file_format="mmCif")

if __name__== "__main__":
    """executes program. If you already have a dataset, you can load it via pickle"""

    # get_data() # used to retrieve sequences b/c advanced PDB search not currently working...

    # DONE : follow the instructions in get_PDB_info()
    # dataset = generate_dataset()

    #Use this line to avoid re-doing the parsing.
    alphads = pickle.load(open(os.getcwd() + "/alpha_ML_ready_dataset.pickle", "rb"))
    betads = pickle.load(open(os.getcwd() + "/beta_ML_ready_dataset.pickle", "rb"))
    aandbds = pickle.load(open(os.getcwd() + "/AandB_ML_ready_dataset.pickle", "rb"))
    adashbds = pickle.load(open(os.getcwd() + "/AdashB_ML_ready_dataset.pickle", "rb"))
    gpcrds = pickle.load(open(os.getcwd() + "/GPCR_ML_ready_dataset.pickle", "rb"))

    # ########################################################################################
    # # TODO : fill in split_dataset()
    predict_intramolecular_contacts(alphads, betads, aandbds, adashbds, gpcrds)
    #######################################################################################
