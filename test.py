
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

list_of_alignments = [["GTT", "GAT"], ["GT-T", "GTAT"], ["GTT-", "GTTC"]]

print(find_multiple_alignment(list_of_alignments))