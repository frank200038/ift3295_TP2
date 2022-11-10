# function from Marc feely
def createMatrix(numRow,numColumn):
    result = [None] * numRow
    for i in range(numRow):
        result[i] = [0]*numColumn
    return result

dic = {"A":"U","U":"A","G":"C","C":"G"}

def get_inverse_complementary_sequence(sequence):
    inverse_sequence = sequence[::-1]
    complementary_sequence = ""
    for i in range(len(inverse_sequence)):
        complementary_sequence += dic[inverse_sequence[i]]
    return complementary_sequence

S = "GCGUGCUUGCGUGCACG"
S_inv_com = get_inverse_complementary_sequence(S)
print(S)
print(S_inv_com)

def alignement_global(A,B):
    lenA = len(A)
    lenB = len(B)
    #initialisation: lenA rows ,lenB cols
    T = createMatrix(lenA + 1, lenB + 1)
    T[0][0] = 0
    for i in range(lenA+1):
        T[i][0] = 0
    for j in range(lenB+1):
        T[0][j] =0  
    #print B
    print("?    ",end = "")
    for i in range(1,lenB+1):
        print(B[i-1],end = "  ")
    print("")
    #firstrow
    print(" "+str(T[0]))
    #fill the table:
    for i in range(1,lenA+1):
        for j in range(1,lenB+1):
            base = max([T[i-1][g] for g in range(j)])

            delta = 1 if (A[i-1] == B[j-1]) else 0
            #print("S:", B[j-1], "T:", A[i-1], "delta:", delta)

            T[i][j] =  base + delta
        #print each row
        print(A[i-1]+str(T[i]))

alignement_global(S_inv_com,S)

#      G  C  G  U  G  C  U  U  G  C  G  U  G  C  A  C  G  
#  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# C[0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0]
# G[0, 1, 0, 2, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2]
# U[0, 0, 1, 1, 3, 2, 2, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2]
# G[0, 1, 0, 2, 1, 4, 3, 3, 3, 4, 3, 4, 3, 4, 3, 3, 3, 4]
# C[0, 0, 2, 1, 2, 2, 5, 4, 4, 4, 5, 4, 4, 4, 5, 4, 5, 4]
# A[0, 0, 0, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 6, 5, 5]
# C[0, 0, 1, 0, 2, 2, 3, 2, 5, 5, 6, 5, 5, 5, 6, 5, 7, 6]
# G[0, 1, 0, 2, 1, 3, 2, 3, 3, 6, 5, 7, 6, 7, 6, 6, 6, 8]
# C[0, 0, 2, 1, 2, 2, 4, 3, 3, 3, 7, 6, 7, 7, 8, 7, 8, 7]
# A[0, 0, 0, 2, 2, 2, 2, 4, 4, 4, 4, 7, 7, 7, 7, 9, 8, 8]
# A[0, 0, 0, 0, 2, 2, 2, 2, 4, 4, 4, 4, 7, 7, 7, 8, 9, 9]
# G[0, 1, 0, 1, 0, 3, 2, 2, 2, 5, 4, 5, 4, 8, 7, 7, 8, 10]
# C[0, 0, 2, 1, 1, 1, 4, 3, 3, 3, 6, 5, 5, 5, 9, 8, 9, 8]
# A[0, 0, 0, 2, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 6, 10, 9, 9]
# C[0, 0, 1, 0, 2, 2, 3, 2, 4, 4, 5, 4, 6, 6, 7, 6, 11, 10]
# G[0, 1, 0, 2, 1, 3, 2, 3, 3, 5, 4, 6, 5, 7, 6, 7, 7, 12]
# C[0, 0, 2, 1, 2, 2, 4, 3, 3, 3, 6, 5, 6, 6, 8, 7, 8, 7]

# \begin{table}[]
# \begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}
# \hline
#   &   & G & C & G & U & G & C & U & U & G & C & G & U & G & C & A & C & G \\ \hline
#   & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \hline
# C & 0 & 0 & 1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 1 & 0 & 1 & 0 \\ \hline
# G & 0 & 1 & 0 & 2 & 1 & 2 & 1 & 1 & 1 & 2 & 1 & 2 & 1 & 2 & 1 & 1 & 1 & 2 \\ \hline
# U & 0 & 0 & 1 & 1 & 3 & 2 & 2 & 3 & 3 & 2 & 2 & 2 & 3 & 2 & 2 & 2 & 2 & 2 \\ \hline
# G & 0 & 1 & 0 & 2 & 1 & 4 & 3 & 3 & 3 & 4 & 3 & 4 & 3 & 4 & 3 & 3 & 3 & 4 \\ \hline
# C & 0 & 0 & 2 & 1 & 2 & 2 & 5 & 4 & 4 & 4 & 5 & 4 & 4 & 4 & 5 & 4 & 5 & 4 \\ \hline
# A & 0 & 0 & 0 & 2 & 2 & 2 & 2 & 5 & 5 & 5 & 5 & 5 & 5 & 5 & 5 & 6 & 5 & 5 \\ \hline
# C & 0 & 0 & 1 & 0 & 2 & 2 & 3 & 2 & 5 & 5 & 6 & 5 & 5 & 5 & 6 & 5 & 7 & 6 \\ \hline
# G & 0 & 1 & 0 & 2 & 1 & 3 & 2 & 3 & 3 & 6 & 5 & 7 & 6 & 7 & 6 & 6 & 6 & 8 \\ \hline
# C & 0 & 0 & 2 & 1 & 2 & 2 & 4 & 3 & 3 & 3 & 7 & 6 & 7 & 7 & 8 & 7 & 8 & 7 \\ \hline
# A & 0 & 0 & 0 & 2 & 2 & 2 & 2 & 4 & 4 & 4 & 4 & 7 & 7 & 7 & 7 & 9 & 8 & 8 \\ \hline
# A & 0 & 0 & 0 & 0 & 2 & 2 & 2 & 2 & 4 & 4 & 4 & 4 & 7 & 7 & 7 & 8 & 9 & 9 \\ \hline
# G & 0 & 1 & 0 & 1 & 0 & 3 & 2 & 2 & 2 & 5 & 4 & 5 & 4 & 8 & 7 & 7 & 8 & 10 \\ \hline
# C & 0 & 0 & 2 & 1 & 1 & 1 & 4 & 3 & 3 & 3 & 6 & 5 & 5 & 5 & 9 & 8 & 9 & 8 \\ \hline
# A & 0 & 0 & 0 & 2 & 2 & 2 & 2 & 4 & 4 & 4 & 4 & 6 & 6 & 6 & 6 & 10 & 9 & 9 \\ \hline
# C & 0 & 0 & 1 & 0 & 2 & 2 & 3 & 2 & 4 & 4 & 5 & 4 & 6 & 6 & 7 & 6 & 11 & 10 \\ \hline
# G & 0 & 1 & 0 & 2 & 1 & 3 & 2 & 3 & 3 & 5 & 4 & 6 & 5 & 7 & 6 & 7 & 7 & 12 \\ \hline
# C & 0 & 0 & 2 & 1 & 2 & 2 & 4 & 3 & 3 & 3 & 6 & 5 & 6 & 6 & 8 & 7 & 8 & 7 \\ \hline
# \end{tabular}
# \end{table}
