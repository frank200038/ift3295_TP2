from sys import argv

'''
Author: YuChen Hui, Yan Zhuang

To run the program, you need to provide the following arguments:

1. The path to the file containing the sequences

Example:
python3 Partie1.py sequence.fq

Here we assume that the fastq file contains only two sequences. (The program will not work if there are more than two 
sequences in the file.)

To change the match, missmatch and indel values, you need to change the values in the main function.
'''

class Cell():
    def __init__(self):
        self.score = 0
        self.provenance = []  # should be a list of an ordered pair (i,j)
        self.position = (0, 0)  # should be an ordered pair (i,j)

    def set_score(self, score):
        self.score = score

    def add_provenance(self, i, j):
        self.provenance.append((i, j))

    def set_position(self, i, j):
        self.position = (i, j)

    def __str__(self) -> str:
        return "[score: {}, provenance: {}, position: {}]".format(self.score, self.provenance, self.position)


def createMatrix(numRow, numColumn):
    result = [None] * numRow
    for i in range(numRow):
        result[i] = []
        for j in range(numColumn):
            cell = Cell()
            cell.set_position(i, j)
            result[i].append(cell)

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

def alignment_prefix_suffix(A, B, match, missmatch, indel, horizontal):
    '''
    If not horizontal In the table of dynamic programming, A is the suquence of column and B is the sequence of row.
    If horizontal In the table of dynamic programming, A is the suquence of row and B is the sequence of column. 
    indel should be negative.

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

    T = createMatrix(lenA + 1, lenB + 1)
    # T[0][0] = Cell()

    if horizontal:
        # if horizontal
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
        # fist column = 0
        # fist row = 0, indel ,2*indel...
        for i in range(lenA + 1):
            pass
            # pass means T[i][0] = Cell()
        for j in range(lenB + 1):
            T[0][j].set_score(indel * j)
            if j != 0:
                # j = 0 -> T[0][0], no provenance
                T[0][j].add_provenance(0, j - 1)

    for i in range(1, lenA + 1):
        for j in range(1, lenB + 1):
            left = T[i][j - 1].score
            diagonal = T[i - 1][j - 1].score
            up = T[i - 1][j].score
            # 
            delta = match if (A[i - 1] == B[j - 1]) else missmatch
            max_score = max(left + indel, up + indel, diagonal + delta)
            T[i][j].set_score(max_score)
            # for tracing back
            if max_score == left + indel:
                T[i][j].add_provenance(i, j - 1)
            if max_score == up + indel:
                T[i][j].add_provenance(i - 1, j)
            if max_score == diagonal + delta:
                T[i][j].add_provenance(i - 1, j - 1)

    ## find the max score in the last row or last column
    if horizontal:
        # find the max score in the last column
        optimal_cell = T[0][lenB]

        for i in range(1, lenA + 1):
            if T[i][lenB].score > optimal_cell.score:
                optimal_cell = T[i][lenB]


    else:
        # find the max score in the last row
        optimal_cell = T[lenA][0]

        for j in range(1, lenB + 1):
            if T[lenA][j].score > optimal_cell.score:
                optimal_cell = T[lenA][j]

    return T, optimal_cell


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


'''
Even though this function will extract all the sequences in a fastq file, we assume that the FASTQ file will contain only
two sequences.
'''

def readFile(path):
    with open(path, 'r') as file:
        lines = file.readlines()

        sequences = []

        count = 0
        for line in lines:
            if count % 4 == 1:
                sequences.append(line)
            count += 1

    return sequences


def main():
    sequences = readFile(argv[1])

    # We assume here we only have to sequences
    # Take first and second sequence as default.
    A = sequences[0].strip()
    B = sequences[1].strip()
    table = alignment_prefix_suffix(A=A, B=B, match=+4, missmatch=-4, indel=-8, horizontal=False)
    # print("The table of dynamic programming is :")
    # print_table_score(table[0])
    print("-----------------------------------------------")
    print("The optimal score is :")
    print(table[1].score)
    print("-----------------------------------------------")
    print("All the optimal paths are:")
    paths = find_optimal_path(table[0], table[1])
    print_paths(paths)
    print("-----------------------------------------------")
    print("All the optimal alignments are:")
    alignments = get_alignment(paths, A, B)
    print_alignments(alignments)
    print("-----------------------------------------------")


if __name__ == "__main__":
    main()
