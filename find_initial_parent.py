#!/usr/bin/env python
''' Author : Huy Nguyen
    Project: provide method and helper method to find all possible parents
             from 2 given gene_block
             Also calculate the cost for each of the possible parents based
             on a given cost function
    Type   : Genome type is a combination of gene type and a '|' character
             Gene type is string stype
    Start  : 01/07/2016
    End    :
'''

#######################################################################################
# Helper functions to count dupplication, split, parse Genome to create
# set of genes, set of blocks of gene
#######################################################################################

''' @function   : return number of split event from a genome
    @para_type  : list (this list comprises of genes, and '|' as split event)
    @result_type: integer
'''
def countSplit(Genome):
    count= 0
    for element in Genome:
        if element == '|':
            count +=1
    return count

''' @function   : return list of genes that are duplicated in a block of gene with the
                  block it appear in
    @para_type  : list (this list comprises of genes, and '|' as split event)
    @result_type: list of gene 
'''
def countDup(Genome):
    # splitting Genome into blocks of gene 
    gene_list=Genome.split('|')
    gene_dup=[]
    # iterate through the list of gene blocks
    for index in range(len(gene_list)):
        # iterate through each gene in a gene block
        for item in range (len(gene_list[index])):
            # iterate from the index next to our gene (need checking) to the end of the gene block
            for i in range(item+1, len(gene_list[index])):
                if gene_list[index][item] == gene_list[index][i]:
                    gene_dup.append(gene_list[index][item])                    
    return gene_dup

''' @function   : return set of gene in a Genome
    @para_type  : list (this list comprises of genes, and '|' as split event)
    @result_type: set of gene
'''
def setOfGene(Genome):
    Genome_gene= set()
    for gene in Genome:
        if gene != '|':
            Genome_gene.add(gene)
    return Genome_gene

''' @function   : return set of geneblocks/gene in a Genome, given that Genome has split event
    @para_type  : list (this list comprises of genes, and '|' as split event)
    @result_type: set of gene blocks( at least exist 1 gene block)
'''
def setOfBlocks(Genome):
    GeneBlocks = Genome.split('|')
    initial = set()
    # sorted before add to make sure 'ab' same as 'ba'
    for element in GeneBlocks:
        initial.add(''.join(sorted(element)))
    return initial

#######################################################################################
# Helper functions to do the reduction
#######################################################################################

''' @function   : return list of blocks genes or genes that has no element that is a subset of an
                  element in the set
    @para_type  : set of block of genes/ genes.
    @result_type: list of blocks genes, genes.
'''
def reductionSubset(initial):
        ''' for cases that return set {'abc','bc','a','ef','f'}
            I will remove 'bc', 'a' and 'f' because having 'abc','ef' is 
            representative enough. In addition, the cost of split is the minimum
            (need to be proved)
        '''
        CandidateBlocks=[]
        # make a copy of our initial set into a list to use index'
        for element in initial:
            CandidateBlocks.append(element)
        # clear our initial set, prepare to add those that are not subsets to it'
        initial.clear()
        for index in range(len(CandidateBlocks)):
            add = True # to add the element or not
            for i in range(len(CandidateBlocks)): # 'ab' might appear before 'b'
                # for the case where the element is only a letter
                if index !=i:
                    if CandidateBlocks[index] in CandidateBlocks[i]:
                        add = False
                        break
                # for cases like 'jk' and 'jhk'
                    subcount = 0
                    for gene in CandidateBlocks[index]:
                        if gene in CandidateBlocks[i]:
                            subcount +=1
                    if subcount == len(CandidateBlocks[index]):
                        add = False
                        break
            if add:
                initial.add(CandidateBlocks[index])
        return initial
    
''' @function   : return list of blocks genes or genes that for each block of genes, the gene
                  in the block has count of 2.
    @para_type  : set of block of genes/ genes.
    @result_type: list of blocks genes/ genes 
'''
def reductionCount(initial, dic):
    result= set()
    for element in initial:
        # create a copy for each element starting as empty
        copy =''
        for gene in element:
            # adding the gene in the element if the count is 2
            if dic[gene][0] == 2:
                copy += gene
        if len(copy) >0:
            result.add(copy)
    return result

#######################################################################################
# Helper functions to take care of transition state of a gene count
# depends on the second parameter (either the target is a genome, or a set)
#######################################################################################
''' @function   : Define the specification for each gene_frequency
    @para_type  : 1 float
    @result_type: 1 integer
'''           
def frequency(gene_frequency):
    if gene_frequency < .25:
        return 0
    elif gene_frequency < .5:
        return 1
    else:
        return 2
''' @function   : Given counts of specific gene g, return the count in their closest common
                  ancestor. this is for set vs Genome
    @para_type  : 2 integers, 1 float
    @result_type: integer
'''
def transitionSG(count1,count2,gene_frequency):
    # for case that only yield 1 answer
    if count1 == 0 and count2 == 0:
        return 0
    elif count1 == 2 and count2 == 1:
        return 2
    else:
        return frequency(gene_frequency)

''' @function   : Given counts of specific gene g, return the count in their closest common
                  ancestor. this is for set vs set
    @para_type  : 2 integers, 1 float
    @result_type: integer
'''
def transitionSS(count1,count2,gene_frequency):
    # for case that only yield 1 answer
    if count1 == 0 and count2 == 0:
        return 0
    elif count1 == 2 and count2 == 2:
        return 2
    else:
        return frequency(gene_frequency)
    
#######################################################################################
# Main functions to find initial set of genes/ genes block
# 3 functions to deal with 3 type of parameter (Genome,Genome), (Set,Genome), (Set,Set)
#######################################################################################
    
''' @function   : find the initial set of blocks of genes/ genes, and provide dictionary that
                  has key is the gene , and value is of type list. The list has 2 index. Index 0
                  is of value 1 or 2 : 1 means it appear in 1 genome, 2 means it appear in both genome.
                  Index 1 is the numer of time the gene appear (this is not different at Genome level,
                  but for set it will be incremented for better correctness)
                  Also return the number of genome it runs through,
                  which is 2
    @para_type  : 2 list of genes , 2 integer number
    @result_type: tuple og (list of gene, dictionary,integer)
'''
def findSetInitial_GG(Genome1,Genome2,split1,split2):
    # set of initial gene and blocks that will be return
    initial =set()
    # element count is a dictionary
    elementCount= {}
    # set of different gene in Genome1
    Genome1_gene= setOfGene(Genome1)
    # set of different gene in Genome2
    Genome2_gene= setOfGene(Genome2)

    # union the above 2 set to get list of possible gene from 2 Genome
    union = Genome1_gene.union(Genome2_gene)
    # intersect the above 2 set to get list of gene that appears in both genome
    intersect = Genome1_gene.intersection(Genome2_gene)
    # create a set of element count
    for gene in union:
        # add a list (gene,value) (since set does not take dictionary)
        if gene in intersect:
            elementCount[gene]=[2,2]
        else:
            elementCount[gene]=[1,1]

    ### if neither of them has a split (life is real good)
    if split1 == 0 and split2 == 0:
        initial = union

    ### if both of them have split ( grrrr )
    if split1 !=0 and split2 !=0:
        ''' parse the genome1 and genome2 into 2 sets og group of genes
            that is delimitered by |
            ex: ab|cd|ef will be same as ef|cd|ab
            ab|cd|ef will be the same as dc|ba|fe
        '''
        # add the set of gene blocks from genome1
        initial.update(setOfBlocks(Genome1))
        # add the set of gene blocks from genome2
        initial.update(setOfBlocks(Genome2))
        # reduce by the count
        initial = reductionCount(initial, elementCount)
        # reduce the subset
        initial = reductionSubset(initial)
        
    # if one of them has a split (we can assume second one has split
    # because we can switch from the parameter of our function)
    if split1 == 0 and split2 !=0:
        # add the set of gene blocks from genome2
        initial.update(setOfBlocks(Genome2))
        # add element that is in genome1 but not 2
        initial=initial.union(Genome1_gene)
        # reduce by the count
        initial = reductionCount(initial, elementCount)
        # reduce the subset
        initial = reductionSubset(initial)
    return (initial,elementCount,2)

''' @function   : find the initial set of blocks of genes/ genes, and provide dictionary that
                  has key is the gene , and value is either 0, 1, 2. 1 means it appears in 1 of them
                  (either my set, or the genome), 2 means it appear in the set and the genome.
                  0 means no appearance.
    @para_type  : 1 list of genes , 1 tuple of (set,dictionary,Genome count), 1 split.
    @result_type: tuple of (set, dictionary,Genome count)
'''
def findSetInitial_SG(myTuple,Genome,split):
    ### working on the Genome info
    # get the gene set from the Genome.
    GenomeDic={}
    Genome_gene= setOfGene(Genome)
    # get the geneblock set from the Genome
    Genome_geneBlocks = setOfBlocks(Genome)
    # create dictionary for the genome Gene
    Genome_dic={}
    for gene in Genome_gene:
        Genome_dic[gene]=1

    ### extract info from myTuple
    # the gene/ gene block set
    initial= myTuple[0]
    # the gene Count dictionary
    elementCount = myTuple[1]
    # increment the size of Leaf(x)
    count= myTuple[2]+1 # increment the count by 1

    ### editting the element count Dic for the new Tuple:

    # case where the Genome.count(g) is 1:
    # create new key for the new gene from Genome, if existing, just update
    for key in Genome_dic:
        if key not in elementCount:
            # this means that for the Set, none of its leave have the gene until this Genome
            elementCount[key]=[0,1]
        else:
            # increment the number of leaf that has this gene (key) by one
            # in elemtCount dictionary
            elementCount[key][1] += 1
            # calculate the frequency
            frequency = elementCount[key][1]/ count
            # use the transitional function helper to update the key value, count of gene in
            # Genome is 1 
            newValue = transitionSG(elementCount[key][0],1,frequency)
            # update the elementCount
            elementCount[key][0]=newValue
    # case where the Genome.count(g) is 0:
    for key in elementCount:
        if key not in Genome_dic:
            # calculate the frequency
            frequency = elementCount[key][1]/ count
            # use the transitional function helper to update the key value, count of gene in
            # Genome is 0
            newValue = transitionSG(elementCount[key][0],0,frequency)
            # update the elementCount
            elementCount[key][0]=newValue

    ### create the initial Set for the new Tuple:
                
    # edit the initialSet from the Set (reducecount)
    initial = reductionCount(initial,elementCount)
    # edit the GenomeBlocks
    Genome_geneBlocks = reductionCount(Genome_geneBlocks,elementCount)
    
    # union the above 2, then do reductionSubset
    resultInitial = initial.union(Genome_geneBlocks)
    resultInitial = reductionSubset(resultInitial)

    return (resultInitial, elementCount, count)
''' @function   : find the initial set of blocks of genes/ genes, and provide dictionary that
                  has key is the gene , and value is either 0, 1, 2. 1 means it appears in 1 of them
                  , 2 means it appear in the set and the genome. 
    @para_type  : 2 tuples of (set,dictionary,Genome count).
    @result_type: tuple of (set, dictionary,Genome count)
'''
def findSetInitial_SS(myTuple1,myTuple2):
    ### extract info from myTuple1    
    # the gene/ gene block set
    initial1= myTuple1[0]
    # the gene Count dictionary
    elementCount1 = myTuple1[1]

    ### extract info from myTuple2
        # the gene/ gene block set
    initial2= myTuple2[0]
    # the gene Count dictionary
    elementCount2 = myTuple2[1]

    ### count of Genome
    count = myTuple1[2] + myTuple2[2]
    ### editting the element count Dic for the new Tuple:

    # create new elementCountdic:
    elementCount ={}
    # case where the initial2.count(g) exist (elemntCount2[key][1] != 0)
    # create new key for the new gene from Genome, if existing, just update
    for key in elementCount2:
        if key not in elementCount1:
            # this means that for the Set, its leaf has none of gene g.
            # calculate frequency:
            frequency = elementCount2[key][1]/count
            # use transition function to derive new value for the count of Set(x)
            newValue = transitionSS(elementCount2[key][0],0,frequency)
            # update the elementCount:
            elementCount[key]=[newValue,elementCount2[key][1]]
        else:
            # increment the number of leaf for the parent node by the sum
            numberOfLeaf = elementCount1[key][1] +elementCount2[key][1]
            # calculate the frequency
            frequency = numberOfLeaf/ count
            # use the transitional function helper to update the key value, count of gene in
            # Genome is 1 
            newValue = transitionSS(elementCount1[key][0],elementCount2[key][0],frequency)
            # update the elementCount
            elementCount[key]=[newValue,numberOfLeaf]
    # case where the Genome.count(g) is 0:
    for key in elementCount1:
        if key not in elementCount2:
            # calculate the frequency
            frequency = elementCount1[key][1]/ count
            # use the transitional function helper to update the key value, count of gene in
            # Genome is 0
            newValue = transitionSS(elementCount1[key][0],0,frequency)
            # update the elementCount
            elementCount[key]=[newValue,elementCount1[key][1]]

    ### create the initial Set for the new Tuple:
                
    # edit the initialSet from the Set (reducecount)
    initial1 = reductionCount(initial1,elementCount)
    initial2 = reductionCount(initial2,elementCount)
    # union both of them
    initial = initial1.union(initial2)
    # edit the GenomeBlocks
    initial = reductionSubset(initial)

    return (initial, elementCount, count)
    


