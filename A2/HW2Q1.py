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
if not sys.warnoptions:
    warnings.simplefilter("ignore")

from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from imblearn.combine import SMOTETomek
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from copy import deepcopy

def calc_residue_dist(residue_one, residue_two) :
    """Computes and returns the distance between two residues, by comparing the position of their alpha carbons"""
    #TODO : return an integer representing the distance between the two residues, in Angstrom
    ca1 = residue_one["CA"]
    ca2 = residue_two["CA"]
    distance = ca1 - ca2
    return abs(distance)

def compute_distance_matrix(residues) :
    """Computes a matrix of size len(Seq) * len(Seq) with distances between each residue."""
    #TODO : return a numpy 2D array of distances between each residue of the structure.
    #Tip : you might want to make sure you filter out consecutive residues at this step.
    ########################################################################################
    matrix = np.zeros((len(residues), len(residues)), np.float)
    for row, residue_one in enumerate(residues):
        for col, residue_two in enumerate(residues):
            ##if two amino acids are within 5 angstroms of each other in 3D, but distant of at least 10 in sequence, the table should have True, else False.
            if abs(row-col) < 10: matrix[row][col] = -1
            else: matrix[row][col] = calc_residue_dist(residue_one, residue_two)

    return matrix

def extract_residues(model):
    """Returns a list of protein residues given a PDB model"""
    #TODO : return a list of protein residues given a PDB model
    residues = []
    for residue in model.get_residues():
        # PDB.is_aa returns true if residue object/string is an amino acid 
        #use standard=True to check 20 AAs (else checks if 3-letter code)
        if (PDB.is_aa(residue, standard=True)):
            residues.append(residue)
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

def run_tmhmm(pdb_file, filename):
    """Runs tmhmm on input fasta files, sends the output to STDout, processes output into secondary structure"""
    #TODO : run TMHMM and output its resulting secondary structure.
    #If you cannot get tmhmm to run, you can write all your sequences to a FASTA file,
    # use the webserver and parse the output file with this function. Otherwise, you can use this function to run
    #TMHMM on some FASTA file and parse its output.
    
    #NOTE: ran tmhmm on all.fasta using webserver
    # in this function, parsing through output file (called tmhmm_results.txt)
    df = pd.read_csv(filename, sep='\t', header=None)
    col = df.loc[df[0] == pdb_file]
    length = int(col[1].reset_index(drop=True)[0][4:])
    idx = col[5].reset_index(drop=True)[0][9:]
    idx = re.split('o|i|-',idx) 
    # NOTE: first position = 1 in TMHMM (NOT 0)
    idx.pop(0)
    idx.pop(-1)
    ss = ''
    prev = 0
    num = 0
    if idx:
        for i in idx:
            num+=1
            if num%2 == 0:
                ss+=('H'*(int(i)-prev+1))
            else:
                ss+=('C'*(int(i)-prev-1))
            prev = int(i)
        # last index until end
        ss+=('C'*(length-int(idx[-1])))
    else:
        ss+=('C'*length)
    return ss

def get_contact_numbers(contact_map):
    """Returns the proportion of residues involved in intramolecular contacts"""
    #TODO : a utility function for the number of contacts in the protein, which you might use to make sure your output makes sense
    return numpy.count_nonzero(contact_map) / float(contact_map.size)

def generate_ML_dataset(sequence,dssp_ss,tm_ss,has_contact,DSSP_vector, TMHMM_vector, oracle):
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
    tm_nine = list(split_by_n(tm_ss, 9))
    dssp_nine = list(split_by_n(dssp_ss, 9))

    for i in enumerate(seq_nine):
        if(len(seq_nine[i[0]]) < 9): continue
        DSSP_vector.append(zip(seq_nine[i[0]], dssp_nine[i[0]]))
        TMHMM_vector.append(zip(seq_nine[i[0]], tm_nine[i[0]]))

    i = 0
    while i<len(sequence):
        if i == 0:
            oracle.append(has_contact[i+4])
        else:
            oracle.append(has_contact[i])
        i = i + 9
    return DSSP_vector, TMHMM_vector, oracle

def split_by_n(seq, n):
    '''Divides a sequence into n subunits'''
    while seq:
        yield seq[:n]
        seq = seq[n:]

def get_PDB_info(dir):
    """Extracts sequence, DSSP secondary structure, TMHMM secondary structure and contact information from PDB files in input directory"""

    #the three vectors you are required to fill.
    DSSP_vector, TMHMM_vector, oracle = [],[],[]

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
        #TODO : compute a distance matrix of size len(sequence)*len(sequence) with the distance between each residue
        matrix = compute_distance_matrix(residues)
        
        #TODO : contact map should be a boolean numpy array of the same size as the distance matrix.
        #if two amino acids are within 5 angstroms of each other in 3D, but distant of at least 10 in sequence, the table should have True, else False.
        # use ">" because consecutive distance in matrix set to -1
        contact_map = (matrix>-1)

        #we can then define the list has_contact from the contactmap array.
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
       
        # get dssp secondary structure (if using biopython)
        for aa in sorted(dssp_info.keys()):
            ss = dssp_info[aa][1]
            # H = G is 3-turn, H = 4-turn, I = 5-turn helix
            if (ss in ['H', 'G', 'I']):
                dssp_ss+="H"
            #everything that isnt a helix is considered a coil...
            else:
                dssp_ss+="C"

        #TODO : write the sequence to a fasta file to call TMHMM with it, or to use the webserver
        filename = write_fasta(sequence,PDB_file)

        # already ran tmhmm on compiled fasta file on webserver -> need to parse file
        #TODO : obtain secondary structure prediction for this FASTA file with TMHMM
        tm_ss = run_tmhmm(filename, os.getcwd()+'/tmhmm_results.txt')

        DSSP_vector, TMHMM_vector, oracle = generate_ML_dataset(sequence,dssp_ss,tm_ss,has_contact,DSSP_vector, TMHMM_vector, oracle)
    return DSSP_vector, TMHMM_vector, oracle

#DONE
def generate_dataset():
    """Runs the PDB parsing utility"""
    DSSP_vector, TMHMM_vector, oracle = get_PDB_info(os.path.join("data","pdb"))
    print("successfully parsed PDB data")
    
    #store a pickle of your results to avoid repeating get_PDB_info
    pickle.dump((DSSP_vector, TMHMM_vector, oracle),open("ML_ready_dataset.pickle","wb"))
    return DSSP_vector, TMHMM_vector, oracle

def split_dataset(X, Y):
    """Splits the dataset into training and testing"""
    # TODO : split the X and Y dataset into a reasonable training set and a test set. Your test set should have have 20% of the datapoints.
    # Tip : look up train_test_split with scikit-learn
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.20)
    return X_train, X_test, Y_train, Y_test

#DONE
def format_simple_dataset(vector,solutions):
    """takes as input a vector of sequence strings and a vector of booleans
    outputs a vector of size 9 vectors of tuples and an array of binary numbers"""

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
    X = np.array(formatted_X)
    formatted_Y = np.zeros(training_size)
    for i in range(training_size):
        formatted_Y[i] = contact_truth.index(solutions[i])
    Y = np.array(formatted_Y)
    return X,Y

#DONE
def format_one_hot_dataset(vector,solutions):
    """Takes as input a vector of 9 sequence, ss tuples, outputs a vector of one-hot vectors of size 198 and a binary vector"""
    AA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T", "X"]
    ss = ["H", "C"]
    contact_truth = [False, True]

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
    X = np.array(formatted_X)
    formatted_Y = np.zeros(training_size)
    for i in range(training_size):
        formatted_Y[i] = contact_truth.index(solutions[i])
    Y = np.array(formatted_Y)
    return X,Y

def run_NN_on_sequence(vector, solutions):
    """Trains a scikit-learn basic 4 layers neural network on a sequence vector, outputs results"""
    #X is a the input vector of ML-ready numpy arrays, Y is the corresponding oracle
    vector_copy = deepcopy(vector)
    X, Y = format_simple_dataset(vector_copy, solutions)    
    #TODO : fill in split_dataset
    # X_train, X_test, Y_train, Y_test = split_dataset(X,Y)
    method = SMOTETomek()
    X_resampled, y_resampled = method.fit_resample(X, Y)
    X_train, X_test, Y_train, Y_test = split_dataset(X_resampled, y_resampled)

    clf = MLPClassifier(solver='sgd', alpha=1e-6, hidden_layer_sizes=(156), activation="logistic", shuffle=True,
                        verbose=False, random_state=1, tol=1e-5, max_iter=350)
    clf.fit(X_train, Y_train)
    Y_pred_train = clf.predict(X_train)
    Y_pred_test = clf.predict(X_test)
    print("score, sequence only")
    print("TRAIN SCORES", classification_report(Y_train, Y_pred_train))
    print("CONFUSION MATRIX\n", confusion_matrix(Y_train, Y_pred_train))
    print("TEST SCORES", classification_report(Y_test, Y_pred_test))
    print("CONFUSION MATRIX\n", confusion_matrix(Y_test, Y_pred_test))

def run_NN_for_ss_type(vector,solutions):
    """Trains a scikit-learn basic 4 layers neural network on a sequence-structure vector, outputs results"""
    #X is a the input vector of ML-ready numpy arrays, Y is the corresponding oracle
    X,Y = format_one_hot_dataset(vector,solutions)
    #TODO : fill in split_dataset
    # X_train, X_test, Y_train, Y_test = split_dataset(X,Y)
    method = SMOTETomek()
    X_resampled, y_resampled = method.fit_resample(X, Y)
    X_train, X_test, Y_train, Y_test = split_dataset(X_resampled, y_resampled)

    clf = MLPClassifier(solver='sgd', alpha=1e-4, hidden_layer_sizes=(156), activation="logistic", tol=1e-5,
                        shuffle=True, verbose=False, random_state=1, max_iter=350)

    clf.fit(X_train,Y_train)
    Y_pred_train = clf.predict(X_train)
    Y_pred_test = clf.predict(X_test)
    print("score, ss and sequence")
    print("TRAIN SCORES", classification_report(Y_train, Y_pred_train))
    print("CONFUSION MATRIX\n", confusion_matrix(Y_train, Y_pred_train))
    print("TEST SCORES", classification_report(Y_test, Y_pred_test))
    print("CONFUSION MATRIX\n", confusion_matrix(Y_test, Y_pred_test))

#DONE
def predict_intramolecular_contacts(dataset):
    """Compares neural network results for DSSP and TMHMM secondary structures"""
    DSSP_vector, TM_vector, solutions = dataset
    print("SEQUENCE ONLY RESULTS")
    #use either vector for sequence only prediction
    run_NN_on_sequence(DSSP_vector,solutions)
    print("DSSP RESULTS")
    run_NN_for_ss_type(DSSP_vector,solutions)
    print("TM RESULTS")
    run_NN_for_ss_type(TM_vector,solutions)

    #TODO : analyze your results!

#DONE
def get_data():
    L = ["1D6G", "1FJR", "1MF6", "1NZS", "1WSO", "1Y7J", "1YM7", "1YTV", "2DCO", "2JX0", "2JX4", "2JX9", "2K3U", "2K9P", "2KI9", "2KOE", "2L60", "2L63", "2LMK", "2LNL", "2LQ4", "2M9O", "2N2F", "2N55", "2NZ1", "2O8Z", "2QKH", "2RGN", "2RH1", "2RRS", "2X57", "2XVT", "2XWT", "2Y00", "2YDV", "2YX8", "3EML", "3H3G", "3I8S", "3N7S", "3N93", "3NKQ", "3ODU", "3PBL", "3RZE",
"3T33", "3UGU", "3UON", "3V2Y", "4DJH", "4DLO", "4DLQ", "4EIY", "4I6O", "4IAR", "4IB4", "4JKV", "4JQI", "4K5Y", "4L6R", "4MBS", "4MQW", "4N6H", "4NUU", "4OR2", "4P39", "4P3A", "4PNI", "4PXZ", "4QIM", "4R7V", "4R7X", "4RMK", "4TND", "4U14", "4U15", "4X1H", "4XES", "4XNV", "4XT1", "4Z35", "4ZUD", "5CXV", "5DHG", "5DSG", "5EE7", "5FTT", "5KVM", "5KZZ", "5NDD", "5NX2", "5O9H", "5OLL", "5OTU", "5T04", "5T1A", "5UEN", "5UNF", "5V57", "5VBL", "5VEW",
"5WB1", "5WIU", "5XJM", "5ZKC", "5ZTY", "6AK3", "6B73", "6BD4", "6BQH", "6C1R", "6CM4", "6CRK", "6D1U", "6D27", "6D32", "6DO1", "6F27", "6F3V", "6F3W", "6F3X", "6F3Y", "6FFI", "6FJ3", "6GPS", "6GPX", "6HLP", "6IDX", "6IIU", "6J20", "6M9T"]

    pdbl = PDB.PDBList(verbose=False)

    for pdb in L:
        path_name = pdbl.retrieve_pdb_file(pdb, pdir='data/pdb/', file_format="mmCif")
    print("DONE")

if __name__== "__main__":
    """executes program. If you already have a dataset, you can load it via pickle"""

    # get_data() # used to retrieve sequences b/c advanced PDB search not currently working...

    # DONE : follow the instructions in get_PDB_info()
    # dataset = generate_dataset()

    #Use this line to avoid re-doing the parsing.
    dataset = pickle.load(open(os.getcwd() + "/ML_ready_dataset.pickle", "rb"))

    ########################################################################################
    # TODO : fill in split_dataset()
    predict_intramolecular_contacts(dataset)
    ########################################################################################
