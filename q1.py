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
