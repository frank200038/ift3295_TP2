from sys import argv
import math
import pickle
import itertools
'''
Author: YuChen Hui, Yan Zhuang

To run the program, you need to execute

> python3 q2.py sequence.fq

To change the match, missmatch and indel values, you need to change the values in the main function.
'''

def read_blosum62():
    file = open("BLOSUM62.txt", "r")
    lines = file.readlines()
    file.close()
    matrix = []
    for i in range(1, len(lines)):
        line = lines[i].split()
        matrix.append(line[1:])
    # to int
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] = int(matrix[i][j])
    return matrix

blosum62_matrix = read_blosum62()
# Symmetric matrix

dictionary = {"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, "H":8, "I":9, "L":10, "K":11, "M":12, "F":13, "P":14, "S":15, "T":16, "W":17, "Y":18, "V":19, "B":20, "Z":21, "X":22, "*":23}

def blosum62(a, b):
    # a and b are amino acids
    return blosum62_matrix[dictionary[a]][dictionary[b]]

sequence1 = "MGEIGFTEKQEALVKESWEILKQDIPKYSLHFFSQILEIAPAAKGLFSFLRDSDEVPHNNPKLKAHAVKVFKMTCETAIQLREEGKVVVADTTLQYLGSIHLKSGVIDPHFEVVKEALLRTLKEGLGEKYNEEVEGAWSQAYDHLALAIKTEMKQEES"

sequence2 = "MEKVPGEMEIERRERSEELSEAERKAVQATWARLYANCEDVGVAILVRFFVNFPSAKQYFSQFKHMEEPLEMERSPQLRKHACRVMGALNTVVENLHDPEKVSSVLSLVGKAHALKHKVEPVYFKILSGVILEVIAEEFANDFPPETQRAWAKLRGLIYSHVTAAYKEVGWVQQVPNATTPPATLPSSGP"

sequence3 = "MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLSHGSAQIKGHGKKVVAALIEAANHIDDIAGTLSKLSDLHAHKLRVDPVNFKLLGQCFLVVVAIHHPAALTPEVHASLDKFLCAVGTVLTAKYR"

sequence4 = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG"

sequence5 = "MERLESELIRQSWRAVSRSPLEHGTVLFSRLFALEPSLLPLFQYNGRQFSSPEDCLSSPEFLDHIRKVMLVIDAAVTNVEDLSSLEEYLATLGRKHRAVGVRLSSFSTVGESLLYMLEKCLGPDFTPATRTAWSQLYGAVVQAMSRGWDGE"

sequences = [sequence1, sequence2, sequence3, sequence4, sequence5]


class Cell():
    def __init__(self):
        self.score = 0
        self.provenance = []  # should be a list of Cells
        self.position = (0, 0)  # should be an ordered pair (i,j)
        belongs_to = None

    def set_score(self, score):
        self.score = score

    def add_provenance(self, i,j):
        self.provenance.append((i,j))

    def set_position(self, i, j):
        self.position = (i, j)
    
    def set_belongs_to(self, belongs_to):
        self.belongs_to = belongs_to

    def __str__(self) -> str:
        return "[score: {}, provenance: {}, position: {}]".format(self.score, self.provenance, self.position)


def createMatrix(numRow, numColumn, belongs_to):
    result = [None] * numRow
    for i in range(numRow):
        result[i] = []
        for j in range(numColumn):
            cell = Cell()
            cell.set_position(i, j)
            result[i].append(cell)
            cell.set_belongs_to(belongs_to)

    return result




def alignment_affine_gap(A, B, h =10, s = 1, horizontal = False):
    '''
    If not horizontal In the table of dynamic programming, A is the suquence of column and B is the sequence of row.
    If horizontal In the table of dynamic programming, A is the suquence of row and B is the sequence of column. 
    indel should be negative.

    h: gap open penalty
    s: gap extension penalty

    '''

    if horizontal:
        t = A
        A = B
        B = t

    lenA = len(A)
    lenB = len(B)

    # initialization: lenA rows ,lenB cols
    # the suquence of column decide the number of rows
    # the suquence of row decide the number of cols


    V = createMatrix(lenA + 1, lenB + 1, "V")
    E = createMatrix(lenA + 1, lenB + 1, 'E')
    F = createMatrix(lenA + 1, lenB + 1, "F")
    G = createMatrix(lenA + 1, lenB + 1, "G")

    # T[0][0] = Cell()

    if horizontal:
        # if horizontal laisse tomber pour TP2, lol.
        # fist row = 0
        # fist column = 0, indel ,2*indel...
        for i in range(lenA + 1):
            T[i][0].set_score(indel * i)
            if i != 0:
                # i = 0 -> T[0][0], no provenance
                T[i][0].add_provenance(i - 1, 0)
        for j in range(lenB + 1):
            pass
            # pass means T[0][j] = Cell()


    else:
        # V (0, 0) = 0
        # V (i, 0) = F(i, 0) = −g(i) = −h − is
        # V (0, j) = E(0, j) = −g(j) = −h − js
        # E(i, 0) = F(0, j) = −∞
        # G(i, 0) = G(0, j) = −∞

        V[0][0].set_score(0)

        for i in range(lenA + 1):
            V[i][0].set_score(-h - s * i)
            if i != 0:
                # i = 0 -> ?[0][0], no provenance
                V[i][0].add_provenance(i - 1, 0)
                F[i][0].add_provenance(i - 1, 0)
                E[i][0].set_score(-math.inf)
                G[i][0].set_score(-math.inf)

        for j in range(lenB + 1):
            V[0][j].set_score(-h - s * j)
            if j != 0:
                # j = 0 -> ?[0][0], no provenance
                V[0][j].add_provenance(0, j - 1)
                E[0][j].add_provenance(0, j - 1)
                F[0][j].set_score(-math.inf)
                G[0][j].set_score(-math.inf)


    # recurrence
    # V(i, j) = max{G(i, j), F(i, j),E(i, j)}
    # G(i, j) = V (i − 1, j − 1) + δ(S[i], T [j])
    # F(i, j) = max{F(i − 1, j) − s, V (i − 1, j) − h − s}
    # E(i, j) = max{E(i, j − 1) − s, V (i, j − 1) − h − s}
    for i in range(1, lenA + 1):
        for j in range(1, lenB + 1):

            G[i][j].set_score(V[i - 1][j - 1].score + blosum62(A[i - 1],B[j - 1]))
            F[i][j].set_score(max(F[i - 1][j].score - s, V[i - 1][j].score - h - s))
            E[i][j].set_score(max(E[i][j - 1].score - s, V[i][j - 1].score - h - s))
            V[i][j].set_score(max(G[i][j].score, F[i][j].score, E[i][j].score))

            # V
            if G[i][j].score == V[i][j].score:
                V[i][j].add_provenance(i - 1, j - 1)
            if F[i][j].score == V[i][j].score:
                V[i][j].add_provenance(i - 1, j)
            if E[i][j].score == V[i][j].score:
                V[i][j].add_provenance(i, j - 1)

    ## find the right bottom cell of V 
    if horizontal:
        # on s'en fout
        # # find the max score in the last column
        # optimal_cell = T[0][lenB]

        # for i in range(1, lenA + 1):
        #     if T[i][lenB].score > optimal_cell.score:
        #         optimal_cell = T[i][lenB]
        print("tnnd, not possible horizontal haobuhao")

    else:
        return V[lenA][lenB], V



def find_optimal_path(T, optimal_cell):
    result_paths = []

    # find all possible paths using backtracking

    def find_paths(cell, path):
        if len(path) == 0:
            path = [cell]

        if len(cell.provenance) == 0:
            result_paths.append(path)
            return result_paths

        for i, j in cell.provenance:
            path_copy = path.copy()
            path_copy.insert(0, T[i][j])
            find_paths(T[i][j], path_copy)

    find_paths(optimal_cell, [])
    return result_paths


def get_alignment(paths, A,B):
    '''
    A is the first sequence (column)
    B is the second sequence (row)
    alignment[0] is the alignment of A
    alignment[1] is the alignment of B
    '''
    alignments = []
    for path in paths:
        previous_position = path[0].position
        alignment = ["", ""]
        for cell in path[1:]:
            current_position = cell.position
            if current_position[0] == previous_position[0] + 1 and current_position[1] == previous_position[1] + 1:
                alignment[0] = alignment[0] + A[current_position[0]-1]
                alignment[1] = alignment[1] + B[current_position[1]-1] 
            if current_position[0] == previous_position[0] + 1 and current_position[1] == previous_position[1]:
                alignment[0] = alignment[0] + A[current_position[0]-1]
                alignment[1] = alignment[1] + "-"
            if current_position[0] == previous_position[0] and current_position[1] == previous_position[1] + 1:
                alignment[0] = alignment[0]+ "-" 
                alignment[1] = alignment[1]+ B[current_position[1]-1] 
            previous_position = current_position
        alignments.append(alignment)
    return alignments


def print_table(T):  # for review of provenance
    numRow = len(T)
    numColumn = len(T[0])
    for i in range(numRow):
        for j in range(numColumn):
            print(T[i][j], end=" ")
        print("")


def print_table_score(T):  # for review of score
    numRow = len(T)
    numColumn = len(T[0])
    for i in range(numRow):
        for j in range(numColumn):
            print(T[i][j].score, end=" ")
        print("")


def print_paths(paths):
    for path in paths:
        print("One of the optimal paths:")
        for cell in path:
            print(cell, end=" ")
        print("")


def print_alignments(alignments):
    for alignment in alignments:
        print("One of the optimal alignments:")
        print(alignment[0])
        print(alignment[1])
        print("the length of alignment is", len(alignment[0]))




def find_central_sequence_index(matrix):
    '''
    index begin with 0 
    '''
    min_index = 0
    min_score_sum =  sum(matrix[0])

    for i in range(1, len(matrix)):
        score_sum = sum(matrix[i]) 
        if score_sum < min_score_sum:
            min_score_sum = score_sum
            min_index = i 
    
    return min_index

########################################    
## The most important function in TP2 ##
########################################    

def find_multiple_alignment(list_of_alignments):
    '''
    the structure of list_of_alignments is
    list_of_alignments = [alignment_1, alignment_2, ...]
    where alignment_i = [sequence_central, sequence_i]
    example:

    '''

    multiple_alignment = [list_of_alignments[0][0], list_of_alignments[0][1]]

    for i in range(1, len(list_of_alignments)):
        multiple_alignment.append("")
        pointer = 0
        while pointer < max(len(list_of_alignments[i][0]),len(multiple_alignment[0])):
            ## special case: the end of the sequence
            ## this situation is not handled by the general case:
            '''
                                    |
                                    v
            newly added pair is GT-T- GT-TC
            and the current multiple alignment is
            GT-T
            GA-T
            GTAT
            GT-T?
            '''
            len1 = len(list_of_alignments[i][0])
            len2 = len(multiple_alignment[0])
            if len1 != len2 and pointer == max(len1, len2)-1:
                if len2 == min(len1, len2):
                    for a in range(len(multiple_alignment) - 1):
                        multiple_alignment[a] += "-" 
                    multiple_alignment[-1] += list_of_alignments[i][1][pointer]
                    
                if len1 == min(len1, len2):
                    multiple_alignment[-1] += "-" 
                return multiple_alignment

            if multiple_alignment[0][pointer] == list_of_alignments[i][0][pointer]:
                multiple_alignment[-1] += list_of_alignments[i][1][pointer]
            else:
                if list_of_alignments[i][0][pointer] == "-":
                    '''
                    if the newly added pair is:
                      |
                      v
                    GT-T
                    GTAT
                    and pointer == 2, then the result is:
                    GTT  ---> GT-T
                    GAT  ---> GA-T
                    GTA  ---> GTAT
                    '''
                    for a in range(len(multiple_alignment) - 1):
                        multiple_alignment[a] = multiple_alignment[a][:pointer] + "-" + multiple_alignment[a][pointer:]
                    multiple_alignment[-1] += list_of_alignments[i][1][pointer]
                    
                elif multiple_alignment[0][pointer] == "-":
                    '''
                    if the newly added pair is:
                    GT-T
                    GTAT
                    the previous multiple alignment is:
                       |
                       v
                    GT--T
                    GTACT
                    GTA?
                    and pointer == 3, then the result is:
                    GT--T
                    GTACT
                    GTA-T
                    '''
                    list_of_alignments[i][0] = list_of_alignments[i][0][:pointer] + "-" + list_of_alignments[i][0][pointer:]
                    list_of_alignments[i][1] = list_of_alignments[i][1][:pointer] + "-" + list_of_alignments[i][1][pointer:]
                    multiple_alignment[-1] += "-"
                
            pointer += 1

    return multiple_alignment

class Status():
    def __init__(self, open, side):
        self.open = open
        self.side = side

    
    
def score_pair_fin_gap(sequence_1, sequence_2, h =10, s = 1 ):
    score = 0
    open = False
    side = None 
    state = Status(open, side) 
    for i in range(len(sequence_1)):
        if sequence_1[i] != "-" and sequence_2[i] != "-":
            Status.open = False 
            score += blosum62(sequence_1[i], sequence_2[i])
        elif sequence_1[i] == "-" and sequence_2[i] != "-":
            if state.open == False or state.open == True and state.side == 2:
                score = score -h - s
                state.open = True
                state.side = 1 
            elif state.open == True and state.side == 1:
                score = score - s
        elif sequence_1[i] != "-" and sequence_2[i] == "-":
            if state.open == False or state.open == True and state.side == 1:
                score = score -h - s
                state.open = True
                state.side = 2 
            elif state.open == True and state.side == 2:
                score = score - s
        else:
            #else both are "-", we will not penalize in this case. 
            # print("error, both are gaps")
            pass
    
    return score



## SP score
def score_SP(multiple_alignment):
    score = 0
    # combinatory of 2 sequences from the multiple alignment
    combinations = list(itertools.combinations(multiple_alignment, 2))
    # print("combinations", combinations)
    for pair in combinations:
        score += score_pair_fin_gap(pair[0], pair[1])
    return score


# consensus sequence
def get_consensus_sequence(multiple_alignment):
    consensus_sequence = ""
    for i in range(len(multiple_alignment[0])):
        column = [sequence[i] for sequence in multiple_alignment]
        # print("column", column)
        max_count = 0
        max_letter = None
        for letter in column:
            count = column.count(letter)
            if count > max_count:
                max_count = count
                max_letter = letter
        consensus_sequence += max_letter
    return consensus_sequence


def main():
    matrix = createMatrix(len(sequences), len(sequences), "shutong")
    matrix_V = createMatrix(len(sequences), len(sequences), "shu")

    print("Please wait while the program is calculating the distance matrix...")
    # t1 = "MGEIG"
    # t2 = "GEMEI"
    # optimal_cell, V,E,F,G = alignment_affine_gap( t1, t2,blosum62_matrix, h =10, s = 1, horizontal =  False)
    # print(optimal_cell)
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            if i != j:
                optimal_cell, V = alignment_affine_gap( sequences[i], sequences[j], h =10, s = 1, horizontal =  False)
                matrix[i][j] = optimal_cell.score
                matrix_V[i][j] = V
            else:
                matrix[i][j] = 0
                matrix_V[i][j] = None

        print(matrix[i])
    
    # with open("matrix.pickle", "wb") as f:
    #     pickle.dump(matrix, f)
    # with open("matrix_V.pickle", "wb") as f:
    #     pickle.dump(matrix_V, f)
    # with open("matrix.pickle", "rb") as f:
    #     matrix = pickle.load(f)
    # with open("matrix_V.pickle", "rb") as f:
    #     matrix_V = pickle.load(f)
    
    # print(matrix)
    

    sequence_central_index = find_central_sequence_index(matrix)
    print("########################################################")
    print("the central sequence is")
    print("########################################################")
    print(sequences[sequence_central_index])
    # find alignment of central sequence and other sequences
    alignment_pairs = []
    for i in range(len(sequences)):
        if i != sequence_central_index:
            V = matrix_V[i][sequence_central_index] 
            optimal_cell = V[len(V) - 1][len(V[0]) - 1]
            paths = find_optimal_path(V, optimal_cell)
            alignments = get_alignment(paths, sequences[i], sequences[sequence_central_index])
            alignment_pairs.append(alignments)
    
    def swap(liste):
        return [liste[1], liste[0]]
    
    # take one alignment from several possible alignments   
    alignment_pairs_simplified = [swap(alignments[0]) for alignments in alignment_pairs]

    print("########################################################")
    print("The 4 alignments are:")
    print("########################################################")
    print(alignment_pairs_simplified)

    multiple_alignment = find_multiple_alignment(alignment_pairs_simplified)
    print("########################################################")
    print("The multiple alignment is:")
    print("########################################################")
    for sequence in multiple_alignment:
        print(sequence)
    
    print("########################################################")
    print("the SP score is", score_SP(multiple_alignment)) 
    print("########################################################")
    print("########################################################")
    print("the consensus sequence is", get_consensus_sequence(multiple_alignment))
    print("########################################################")


    # these sequences are alignment obtained by online bio informatics tool 
    s1 = "----------------MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLS-H------GSAQIKGHGKKVVAALIE------AANHIDDIAGTLSKLSDLHAHKLRVDPVNFKLLGQCFLVVVAIHHPAALTPEVHASLDKFLCAVGTVLTAKYR-----------------------" 
    s2 = "MEKVPGEMEIERRERSEELSEAERKAVQATWARLYANCEDVGVAILVRFFVNFPSAKQYFSQFKHM-EEPLEMERSPQLRKHACRVMGALNTVV---ENLHDPEKVSSVLSLVGKAHALKHKVEPVYFKILSGVILEVIAEEFANDFPPETQRAWAKLRGLIYSHVTAAYKEVGWVQQVPNATTPPATLPSSGP"
    s3 = "----------------MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHL-KSEDEMKASEDLKKHGATVLTALGGIL---KKKGHH---EAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG-----------------"
    s4 = "-------------MGEIGFTEKQEALVKESWEILKQDIPKYSLHFFSQILEIAPAAKGLFSFLRDSDEVPHN---NPKLKAHAVKVFKMTCETAIQLREEGKVVVADTTLQYLGSIHLK-SGVIDPHFEVVKEALLRTLKEGLGEKYNEEVEGAWSQAYDHLALAIKTE-----MKQEES--------------"
    s5 = "------------------MERLESELIRQSWRAVSRSPLEHGTVLFSRLFALEPSLLPLFQYNGRQFSSPEDCLSSPEFLDHIRKVMLVID-AA--VTNVEDLSSLEEYLATLGRKHRA-VGVRLSSFSTVGESLLYMLEKCLGPDFTPATRTAWSQLYGAVVQAMSRG-----WDGE----------------"
    alignment_outil = [s1,s2,s3,s4,s5]
    print("########################################################")
    print("the SP score of the alignment got by online tools is", score_SP(alignment_outil)) 
    print("########################################################")





if __name__ == "__main__":



    main()







