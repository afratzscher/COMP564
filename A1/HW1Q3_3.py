###NOTE:: FIGURE OUT ERROR! -> MSE or truth **2 ?
import sys, re
from collections import Counter
import math
import os
from tabulate import tabulate
from statistics import mean, stdev
# ============================================ Student info methods================================================
def get_student_name():
    # @TO_STUDENT: Write your name here
    student_name = "Anne-Sophie Fratzscher"
    if not student_name:
        raise ValueError("Error: you forgot to add your name in get_student_name method.")
    return student_name

def get_student_id():
    # @TO_STUDENT: Write your student id here
    student_id = "260705446"
    if not student_id:
        raise ValueError("Error: you forgot to add your student id in get_student_id method.")
    return student_id
# =================================================================================================================

# =================================== Validate input/output methods================================================
def validate_Q3_1_input_format(subopt_result):
    if not isinstance(subopt_result, list) or [sr for sr in subopt_result if re.search("[^\(\)\.]", sr)]:
        raise ValueError("Error: your input should be list of strings (secondary structures, ie alphabet is '().')")

def validate_Q3_1_output_format(result):
    if not isinstance(result, list) or [sr for sr in result if not isinstance(sr, list)]:
        raise ValueError("Error: your output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ].")
    if [sr for sr in result if sr[0] >= sr[1]]:
        raise ValueError("Error: check your i and j values. They should satisfy: i > j.")
# =================================================================================================================
# ================================================== Helper methods================================================
def parse_subopt_result_file(filepath):
    '''
    Parsing of a standard txt file that contains result of
        "RNAsubopt -p __k__ < myFasta.fasta > subopt_result_file.txt"
        where __k__ is parameter. (Filename chosen randomly. Please, feel free to use your own names.)
    @args:
    filepath: (full or relative) path to the subopt_result_file.txt.
    (Filename chosen randomly. Please, feel free to use your own names.)

    @return: subopt_result: list of the strings (assumed to be secondary structures)
    '''
    subopt_result = []
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if i < 1:
                continue
            subopt_result.append(line.strip())
    return subopt_result

def parse_dot_ps_file(filepath):
    '''
    Parsing of a dot.ps file that contains result of RNAfold program
    @args:
    filepath: (full or relative) path to the dot.ps.
    @return:
    dot_ps_result: list f lists with i, j, freq_i_j
    '''
    dot_ps_result = []
    with open(filepath, 'r') as f:
        is_data = False
        for line in f:
            if not is_data and line.startswith('%start of base pair probability data'):
                is_data = True
                continue
            elif is_data and line.startswith('showpage'):
                break
            elif is_data:
                if line.find('ubox') > 0:
                    # take only first 3 numbers
                    data_line = line.split()[:3]
                    dot_ps_result.append(
                        [int(data_line[0]), int(data_line[1]), float(data_line[2])]
                    )
    return dot_ps_result

# =================================================================================================================
def get_answer_Q3_1(subopt_result):
    '''
    This method should be implemented by student.
    @args:
    subopt_result: a list of the secondary structures produced by RNAsubopt -p <k> for particular input

    @return:
    result: list of lists (example is below) with indices and relevant frequency.
    example [ [0, 1, 0.10], [0, 2, 0.15], [0, 3, 0.16], ....... ]

    @note:
    check input/output as advised in code. Question will be marked as 0 in case of not following the formats.
    '''
    # basic check for the proper input
    validate_Q3_1_input_format(subopt_result)
    # @TO_STUDENT: Write your code here
    # @TO_STUDENT: use result variable for results. below is an example of an expected format for result.
   
    result = []
    stack = []
    pairs = []    

    # file already read by parse_subopt_result_file -> each secondary struct. given as string in subopt_resultsc
    for substr in subopt_result:
        for idx, x in enumerate(substr): # idx = index, x = '(', '.', or ')'
        # idx starts at 0 (for first)
            if x == '(':
                stack.append(idx)
            elif x == ')':
                temp = stack.pop()
                pairs.append(tuple([temp+1,idx+1])) # store as tuple so that can use Counter 
                # index adjusted so that 1st nucleotide = idx 1
            else:
                continue

    # count number of occurences using Counter
    counts = Counter(pairs)
    for idx in counts:
        result.append([idx[0], idx[1], float(counts[idx]/len(subopt_result))])
    result = sorted(result, key = lambda x: (x[0], x[1])) # sorts by first index, then second index (e.g. [1,2], [1,4], [2,3] ...)

    # @TO_STUDENT: output should be [ [i1, j1, freq_i1_j1 ], [i2, j2, freq_i2_j2 ], ...  ]
    # use validate_Q3_output_format(result) to validate the output
    validate_Q3_1_output_format(result)
    return result

def get_answer_Q3_2(q3_1_result, dot_ps_result):
    '''
    This method should be implemented by student.
    Compare output from RNAfold and result of question3_1 for the same sequence and return an error (see text assignment)
    result_error is expected to be numeric
    '''
    result_error = 0
    # @TO_STUDENT: Write your code here (trust me, answer is not 0 :-) )

    # freq_tuples = [(x[0], x[1]) for x in q3_1_result]
    # prob_tuples = [(x[0], x[1]) for x in dot_ps_result]
    # # print(freq_tuples)
    # # print(prob_tuples)

    # only_in_prob = (set(prob_tuples) - set(freq_tuples))
    # only_in_freq = (set(freq_tuples) - set(prob_tuples))
    # in_both = (set(prob_tuples) & set(freq_tuples))
    # print(len(only_in_prob))
    # print(len(only_in_freq))
    # print(len(in_both))

    n=0
    for prob in dot_ps_result: # dot_ps_result contains [[i,j, freq]...]
        for freq in q3_1_result:
            if ((freq[0] < freq[1]) and (prob[0] < prob[1])):
                if ((freq[0] == prob[0]) and (freq[1] == prob[1])):
                    n+=1
                    result_error += ((prob[2] **2) - freq[2]) ** 2
                        # b/c values in dot.ps are SQUARE ROOT (so need to square)
    result_error = (math.sqrt(result_error))

    return result_error

# @TO_STUDENT: You can test your methods below by calling methods. Workflow is given already (can be changed).
# @TO_STUDENT: Everything below this point will not be considered as a solution and will be deleted for grading.
# @TO_STUDENT: Advise: be careful with idents. Use only tabs, or only FOUR spaces. NEVER mix them.

print("This is a solution of %s, student_id is %s" % (get_student_name(), get_student_id()) )

run = 1
val_10 = []
val_100 = []
val_1000 = []
val_10000 = []
k_values = [10, 100, 1000, 10000]

while (run <= 10):
    for k in k_values:
        subopt_result_filepath = "subopt" + str(k) + ".txt"
        dot_ps_filepath = "dot.ps"
        command = 'RNAsubopt -p ' + str(k) + ' < HW1Q3.fasta > '+ subopt_result_filepath
        os.system(command)

        subopt_result = parse_subopt_result_file(subopt_result_filepath)
        q3_1_result = get_answer_Q3_1(subopt_result)
        dot_ps_result = parse_dot_ps_file(dot_ps_filepath)
        q3_2_result = get_answer_Q3_2(q3_1_result, dot_ps_result)
        if k == 10:
            val_10.append(q3_2_result)
        elif k == 100:
            val_100.append(q3_2_result)
        elif k == 1000:
            val_1000.append(q3_2_result)
        elif k == 10000:
            val_10000.append(q3_2_result)
        os.remove(subopt_result_filepath)
    run+=1

print('Run             k=10                     k=100                         k=1000                    k=10000')
run = 1
while (run <= 10):
    print(str(run) + '       ' + str(val_10[run-1]) + '         ' + str(val_100[run-1]) + '          ' + str(val_1000[run-1]) + '         ' +       str(val_10000[run-1]))
    run+=1

print('mean     ' + str(mean(val_10)) + '       ' + str(mean(val_100)) + '       ' + str(mean(val_1000)) + '        ' + str(mean(val_10000)))
print('stdev     ' + str(stdev(val_10)) + '       ' + str(stdev(val_100)) + '       ' + str(stdev(val_1000)) + '        ' + str(stdev(val_10000)))
