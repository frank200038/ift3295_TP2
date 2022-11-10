from sys import argv
import math
'''
Author: YuChen Hui, Yan Zhuang

To run the program, you need to execute

> python3 q2.py sequence.fq

To change the match, missmatch and indel values, you need to change the values in the main function.
'''

# Symmetric matrix
blosum62_matrix= [
   [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
   [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4],
   [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4],
   [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4],
   [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4],
   [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4],
   [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
   [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4],
   [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4],
   [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4],
   [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4],
   [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4],
   [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4],
   [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4],
   [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4],
   [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4],
   [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4],
   [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4],
   [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4],
   [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4],
   [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -3, -4, -3, -3, 4, 1, -1, -4],
   [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
   [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4],
   [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]
]

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

    def add_provenance(self, cell):
        self.provenance.append(cell)

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
                alignment[0] = alignment[0] + "-"
                alignment[1] = alignment[1] + B[current_position[1]-1]
            if current_position[0] == previous_position[0] and current_position[1] == previous_position[1] + 1:
                alignment[0] = alignment[0]+ A[current_position[0]-1]
                alignment[1] = alignment[1]+ "-"
            previous_position = current_position
        alignments.append(alignment)
    return alignments


def alignment_affine_gap(A, B, blosum62_matrix, h =10, s = 1, horizontal = False):
    '''
    If not horizontal In the table of dynamic programming, A is the suquence of column and B is the sequence of row.
    If horizontal In the table of dynamic programming, A is the suquence of row and B is the sequence of column. 
    indel should be negative.

    h: gap open penalty
    s: gap extension penalty

    Decide if it is prefix-suffix alignment (horizontal = true) or suffix-prefix alignment(horizontal = false)
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
                V[i][0].add_provenance(V[i - 1][0])
                F[i][0].add_provenance(F[i - 1][0])
                E[i][0].set_score(-math.inf)
                G[i][0].set_score(-math.inf)

        for j in range(lenB + 1):
            V[0][j].set_score(-h - s * j)
            if j != 0:
                # j = 0 -> ?[0][0], no provenance
                V[0][j].add_provenance(V[0][j - 1])
                E[0][j].add_provenance(E[0][j - 1])
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

            if G[i][j].score == V[i][j].score:
                V[i][j].add_provenance(G[i][j])
            if F[i][j].score == V[i][j].score:
                V[i][j].add_provenance(F[i][j])
            if E[i][j].score == V[i][j].score:
                V[i][j].add_provenance(E[i][j])
            
            if F[i-1][j].score -s == F[i][j].score:
                F[i][j].add_provenance(F[i-1][j])
            if V[i-1][j].score - h - s == F[i][j].score:
                F[i][j].add_provenance(V[i-1][j])

            if E[i][j-1].score -s == E[i][j].score:
                E[i][j].add_provenance(E[i][j-1])
            if V[i][j-1].score - h - s == E[i][j].score:
                E[i][j].add_provenance(V[i][j-1])

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
        return V[lenA][lenB], V, E, F, G







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



def main():
    for i in range(len(sequences)):
        matrix = createMatrix(len(sequences), len(sequences), "shutong")

    print("Please wait while the program is calculating the distance matrix...")

    for i in range(len(sequences)):
        for j in range(len(sequences)):
            if i != j:
                optimal_cell, V,E,F,G = alignment_affine_gap( sequences[i], sequences[j],blosum62_matrix, h =10, s = 1, horizontal =  False)
                matrix[i][j] = optimal_cell.score
            else:
                matrix[i][j] = 0

        print(matrix[i])



if __name__ == "__main__":
    main()







