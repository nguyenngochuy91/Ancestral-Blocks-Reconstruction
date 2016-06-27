#!/usr/bin/env python
from __future__ import division
''' Author : Huy Nguyen
    Project: provide method and helper method to find minaml possible parents
             from 2 given gene_block. As well as provide he edit distance
    Type   : Genome type is a combination of gene type and a '|' character
             Gene type is string stype
    Start  : 01/07/2016
    End    : 03/07/2016
'''

#######################################################################################
# Helper functions to count dupplication, split, parse Genome to create
# set of genes, set of blocks of gene
#######################################################################################

''' @function   : return number of split event from a genome
    @input      : list (this list comprises of genes, and '|' as split event)
    @output     : integer
'''
def countSplit(Genome):
    count= 0
    for element in Genome:
        if element == '|':
            count +=1
    return count

''' @function   : return list of genes that are duplicated in a block of gene with the
                  block it appear in
    @input      : list (this list comprises of genes, and '|' as split event)
    @output     : list of gene 
'''
def countDup(Genome):
    # splitting Genome into blocks of gene 
    gene_list=Genome.split('|')
    gene_dup=set()
    # iterate through the list of gene blocks
    for item in gene_list:
        # iterate through each gene in a gene block
        for index in range (len(item)):
            # iterate from the index next to our gene (need checking) to the end of the gene block
            for i in range(index+1, len(item)):
                if item[index] == item[i]:
                    gene_dup.add(item[index])
    return gene_dup

''' @function   : return set of gene in a Genome
    @input      : list (this list comprises of genes, and '|' as split event)
    @output     : set of gene
'''
def setOfGene(Genome):
    Genome_gene= set()
    for gene in Genome:
        if gene != '|':
            Genome_gene.add(gene)
    return Genome_gene

''' @function   : return set of geneblocks/gene in a Genome, given that Genome has split event
    @input      : list (this list comprises of genes, and '|' as split event)
    @output     : set of gene blocks( at least exist 1 gene block)
'''
def setOfBlocks(Genome):
    GeneBlocks = Genome.split('|')
    initial = set()
    # sorted before add to make sure 'ab' same as 'ba'
    for element in GeneBlocks:
        initial.add(''.join(sorted(element)))
    return initial

#######################################################################################
# Helper functions to calculate the edit distance
#######################################################################################
def del_distance(initial1,initial2,initial):
    count =0
    
    initial_gene = set()
    for block in initial:
        gene_set= setOfGene(block)
        initial_gene= initial_gene.union(gene_set)
    initial1_gene = set()
    
    for block in initial1:
        gene_set= setOfGene(block)
        initial1_gene= initial1_gene.union(gene_set)
    initial2_gene = set()
    
    for block in initial2:
        gene_set= setOfGene(block)
        initial2_gene= initial2_gene.union(gene_set)
        
    count += (len(initial_gene-initial1_gene)+len(initial1_gene-initial_gene)+
    len(initial_gene-initial2_gene)+len(initial2_gene-initial_gene))
    return count

def dup_distance(initial1,initial2,initial):
    count =0
    
    initial_gene = set()
    for block in initial:
        gene_set= countDup(block)
        initial_gene= initial_gene.union(gene_set)
        
    initial1_gene = set()   
    if type(initial1) is str: # if it comes from leaf node
        for block in setOfBlocks(initial1):
            dup_gene = countDup(block)
            initial1_gene= initial1_gene.union(dup_gene)
    else: # if it is already in blocks
        for block in initial1:
            dup_gene = countDup(block)
            initial1_gene= initial1_gene.union(dup_gene)
        
    initial2_gene = set()   
    if type(initial2) is str: # if it comes from leaf node
        for block in setOfBlocks(initial2):
            dup_gene = countDup(block)
            initial2_gene= initial2_gene.union(dup_gene)
    else: # if it is already in blocks
        for block in initial2:
            dup_gene = countDup(block)
            initial2_gene= initial2_gene.union(dup_gene)
        
    count += (len(initial_gene-initial1_gene)+len(initial1_gene-initial_gene)+
    len(initial_gene-initial2_gene)+len(initial2_gene-initial_gene))
    return count
    
def split_distance(initial1,initial2,initial):  
    if len(initial1) == 0:
        count= abs(len(initial)-len(initial2))
    else:
        if len(initial2) == 0:
            count= abs(len(initial)-len(initial1))
        else:
            count= abs(len(initial)-len(initial1))+abs(len(initial)-len(initial2))
            
    return count
#######################################################################################
# Helper functions to do the reduction
#######################################################################################

''' @function   : return list of blocks genes or genes that has no element that is a subset of an
                  element in the set
    @input      : set of block of genes/ genes.
    @output     : list of blocks genes, genes.
'''
def reductionSubset(initial,elementCount):
        ''' for cases that return set {'abc','bc','a','ef','f'}
            I will remove 'bc', 'a' and 'f' because having 'abc','ef' is 
            representative enough. In addition, the cost of split is the minimum
            (need to be proved)
        '''
        # findout all the letter that appear in the 
        dic={}
        for item in initial:
            dic[item]={}
            dic[item]['add']=True
            for alphabet in elementCount:
                dic[item][alphabet]=0
            for letter in item:
                dic[item][letter]+=1
        # check for subset

        for item in dic:
            for candidate in dic:
                flag =True 
                if candidate != item:
                    for key in dic[candidate]:
                        if key == "add":
                            continue
                        if dic[item][key]  > dic[candidate][key]:
                            flag = False
           
                            break
                        
                    if flag: # means that all count of letter in item less or 
                            # equal to the candidate => subset
                        dic[item]['add']=False
        initial.clear()
        for item in dic:
            if dic[item]['add']:
                initial.add(item)
        return initial
    
''' @function   : return list of blocks genes or genes that for each block of genes, the gene
                  in the block has count of 2.
    @input      : set of block of genes/ genes.
    @output     : set of blocks genes/ genes 
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
    
''' @function   : return list of blocks genes or genes that for each block of genes, 
                  and depends on the duplication dic to indicate whether to have the
                  or not
    @input      : set of block of genes/ genes.
    @output     : set of blocks genes/ genes 
'''
def reductionDup(initial,duplication):
    result= set()
    for element in initial:
        # create a copy for each element starting as empty
        copy =''
        for index in range(len(element)):
            if index ==0:
                copy += element[index]
            else:
                if element[index] == element[index-1]: #check if the next is dup:
                    if duplication[element[index]][0]==2:
                        copy += element[index]
                else:
                    copy += element[index]
        if len(copy) >0:
            result.add(copy)
    return result


#######################################################################################
# Helper functions to take care of transition state of a gene count
# depends on the second parameter (either the target is a genome, or a set)
#######################################################################################
''' @function   : Define the specification for each gene_frequency
    @input      : 1 float
    @output     : 1 integer
'''           
def frequency(gene_frequency):
    if gene_frequency < .25:
        return 0
    elif gene_frequency <= .5:
        return 1
    else:
        return 2
''' @function   : Given counts of specific gene g, return the count in their closest common
                  ancestor. this is for set vs Genome
    @input      : 2 integers, 1 float
    @output     : integer
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
    @input      : 2 integers, 1 float
    @output     : integer
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
# Helper functions to update element dictionary and duplication dictionary
# depends on the second parameter (either the target is a genome, or a set)
#######################################################################################
        
''' @function   : Given dictionary from a tuple, and a genome. Return the new
                   update dictionary
    @input      : dictionary, dictionary/set
    @output     : dictionary
'''   
def update_dictionary_SG(tuple_dic,genome_dic,count):
    for key in genome_dic:
        if key not in tuple_dic:
            # this means that for the Set, none of its leave have the gene until this Genome
            tuple_dic[key]=[0,1]
        else:
            # increment the number of leaf that has this gene (key) by one
            # in elemtCount dictionary
            tuple_dic[key][1] += 1
            # calculate the frequency
            frequency = tuple_dic[key][1]/ count
            # use the transitional function helper to update the key value, count of gene in
            # Genome is 1 
            newValue = transitionSG(tuple_dic[key][0],1,frequency)
            # update the elementCountfrom __future__ import division
            tuple_dic[key][0]=newValue
    # case where the Genome.count(g) is 0:
    for key in tuple_dic:
        if key not in genome_dic:
            # calculate the frequency
            frequency = tuple_dic[key][1]/ count
            # use the transitional function helper to update the key value, count of gene in
            # Genome is 0
            newValue = transitionSG(tuple_dic[key][0],0,frequency)
            # update the elementCount
            tuple_dic[key][0]=newValue
    return (tuple_dic)
    
''' @function   : Given dictionary from a tuple, and another tuple. Return the new
                   update dictionary
    @input      : dictionary, dictionary/set
    @output     : dictionary
'''   
def update_dictionary_GG(tuple1_dic,tuple2_dic,count):
    result_dic={}
    for key in tuple2_dic:
        if key not in tuple1_dic:
            # this means that for the Set, its leaf has none of gene g.
            # calculate frequency:
            frequency = tuple2_dic[key][1]/count
            # use transition function to derive new value for the count of Set(x)
            newValue = transitionSS(tuple2_dic[key][0],0,frequency)
            # update the elementCount:
            result_dic[key]=[newValue,tuple2_dic[key][1]]
        else:
            # increment the number of leaf for the parent node by the sum
            numberOfLeaf = tuple1_dic[key][1] +tuple2_dic[key][1]
            # calculate the frequency
            frequency = numberOfLeaf/ count
            # use the transitional function helper to update the key value, count of gene in
            # Genome is 1 
            newValue = transitionSS(tuple1_dic[key][0],tuple2_dic[key][0],frequency)
            # update the elementCount
            result_dic[key]=[newValue,numberOfLeaf]
    # case where the Genome.count(g) is 0:
    for key in tuple1_dic:
        if key not in tuple2_dic:
            # calculate the frequency
            frequency = tuple1_dic[key][1]/ count
            # use the transitional function helper to update the key value, count of gene in
            # Genome is 0
            newValue = transitionSS(tuple1_dic[key][0],0,frequency)
            # update the elementCount
            result_dic[key]=[newValue,tuple1_dic[key][1]]
    return result_dic
#######################################################################################
# Functions to find initial set of genes/ genes block
# 3 functions to deal with 3 type of parameter (Genome,Genome), (Set,Genome), (Set,Set)
#######################################################################################
    
''' @function   : find the initial set of blocks of genes/ genes, and provide dictionary that
                  has key is the gene , and value is of type list. The list has 2 index. Index 0
                  is of value 1 or 2 : 1 means it appear in 1 genome, 2 means it appear in both genome.
                  Index 1 is the numer of time the gene appear (this is not different at Genome level,
                  but for set it will be incremented for better correctness)
                  Also return the number of genome it runs through,
                  which is 2
    @input      : 2 list of genes , 2 integer number
    @output     : tuple og (list of gene, dictionary,integer,dictionary)
'''
def findSetInitial_GG(Genome1,Genome2,split1,split2):
    # set of initial gene and blocks that will be return
    initial =set()
    # element count is a dictionary
    elementCount= {}
    # duplication is a dictionary 
    duplication ={}
    # set of different gene in Genome1
    Genome1_gene= setOfGene(Genome1)
    #print('Genome1_gene',Genome1_gene)
    # set of different gene in Genome2
    Genome2_gene= setOfGene(Genome2)
    #blocks of Genome1,Genome2
    Genome1_block= setOfBlocks(Genome1)
    Genome2_block= setOfBlocks(Genome2)

    # union the above 2 set to get list of possible gene from 2 Genome
    union = Genome1_gene.union(Genome2_gene)
    # intersect the above 2 set to get list of gene that appears in both genome
    intersect = Genome1_gene.intersection(Genome2_gene)
    #print('intersect',intersect)
    # create a set of element count
    for gene in union:
        # add a list (gene,value) (since set does not take dictionary)
        if gene in intersect:
            elementCount[gene]=[2,2]
        else:
            elementCount[gene]=[1,1]
    # dealing with duplication
    # for case there is duplication      
    Genome1_dup_dic = countDup(Genome1)
    Genome2_dup_dic = countDup(Genome2)
    union_dup = Genome1_dup_dic.union(Genome2_dup_dic)
    intersect_dup = Genome1_dup_dic.intersection(Genome2_dup_dic)
    # pulling info into the duplication dic
    for gene in union_dup:
        if gene in intersect_dup:
            duplication[gene]=[2,2]
        else:
            duplication[gene]=[1,1]
    ### if neither of them has a split (life is real good)
    if split1 == 0 and split2 == 0:

        string =''
        mylist= list(intersect)
        mylist.sort()
        # for case there is duplication 
        for gene in mylist:
            string+=gene
            if gene in duplication: # check if gene is duplicated at least one)
                if duplication[gene][0]==2: # only if dup is in 2 genome
                    string +=gene
        initial.add(string)
        
    ### if both of them have split ( grrrr )
    if split1 !=0 and split2 !=0:
        ''' parse the genome1 and genome2 into 2 sets og group of genes
            that is delimitered by |
            ex: ab|cd|ef will be same as ef|cd|ab
            ab|cd|ef will be the same as dc|ba|fe
        '''
        # add the set of gene blocks from genome1
        initial.update(Genome1_block)
        # add the set of gene blocks from genome2
        initial.update(Genome2_block)
        
    # if one of them has a split (we can assume second one has splitenome count
    # because we can switch from the parameter of our function)split_total
    if split1 == 0 and split2 !=0:
        # add the set of gene blocks from genome2
        initial.update(Genome2_block)
        # print initial
        new_set= set()
        new_set.add(''.join(sorted(list(Genome1))))
        # add element that is in genome1 but not 2
        initial=initial.union(new_set)
        # print initial

    initial = reductionCount(initial, elementCount)
    if len(duplication)!=0: # only do the reduction Dup if exist duplication 
        initial = reductionDup(initial,duplication)

    initial = reductionSubset(initial,elementCount)

    # find the local cost, and accumulating cost by storing in 2 index list.
    # deletion    
    deletion = del_distance(Genome1,Genome2,initial)    
    deletion_cost=[deletion,deletion]
    # duplication
    if len(duplication)!=0:
        dup = dup_distance(Genome1,Genome2,initial)
    else:
        dup = 0
    duplication_cost=[dup,dup]
    #split
    # reducebyCount for Genome1_blcok
    Genome1_block= reductionCount(Genome1_block,elementCount)
    Genome2_block= reductionCount(Genome2_block,elementCount)
    split = split_distance(Genome1_block,Genome2_block,initial)
    split_cost= [split,split]
    
    return (initial,elementCount,2,duplication,deletion_cost,duplication_cost,split_cost)
    
''' @function   : find the initial set of blocks of genes/ genes, and provide dictionary that
                  has key is the gene , and value is either 0, 1, 2. 1 means it appears in 1 of them
                  (either my set, or the genome), 2 means it appear in the set and the genome.
                  0 means no appearance.
    @input      : 1 list of genes , 1 tuple of (set,dictionary,Genome count), 1 split.
    @output     : tuple of (set, dictionary,integer,dictionary))
'''
def findSetInitial_SG(myTuple,Genome,split):
    ### working on the Genome info
    # get the gene set from the Genome
    Genome_gene= setOfGene(Genome)
    
    # get the geneblock set from the Genome
    Genome_geneBlocks = setOfBlocks(Genome)
    #print('Genome_geneBlocks',Genome_geneBlocks)
    # create dictionary for the genome Gene
    Genome_dic={}
    for gene in Genome_gene:
        Genome_dic[gene]=1

    ### extract info from myTuple
    # the gene/ gene block set
    initial1= myTuple[0]
    #print('initial',initial)
    # the gene Count dictionary
    # increment the size of Leaf(x)
    count= myTuple[2]+1 # increment the count by 1

    ### editting the element count Dic for the new Tuple:

    elementCount= update_dictionary_SG(myTuple[1],Genome_dic,count)
    
    ### create the initial Set for the new Tuple:
    
    # edit the initialSet from the Set (reducecount)
    initial_1 = reductionCount(initial1,elementCount)
    # print('initial',initial)
    # edit the GenomeBlocks
    Genome1_block = reductionCount(Genome_geneBlocks,elementCount)
    # print('Genome_geneBlocks',Genome_geneBlocks)
    # union the above 2, then do reductionSubset
    initial = initial_1.union(Genome1_block)

    ### deal with duplication
    # pull duplication dic from the tuple
    duplication = myTuple[3]
    Genome_dup_dic = countDup(Genome)
    ### editting the duplication Dic for the new Tuple:
    
    duplication= update_dictionary_SG(duplication,Genome_dup_dic,count)
    
    if len(duplication)!=0: # only do the reduction Dup if exist duplication 
        initial = reductionDup(initial,duplication)
    ### final step, reduce subset
    initial = reductionSubset(initial,elementCount)
    
    ### calculate edit distances 
    #deletion
    deletion = del_distance(initial1,Genome,initial)
    accumulate_del = myTuple[4][1] # pull the accumulation del cost
    deletion_cost =[deletion,accumulate_del+deletion]
    #duplication
    if len(duplication)!=0:
        dup = dup_distance(initial1,Genome_geneBlocks,initial)
    else:
        dup =0
    accumulate_dup = myTuple[5][1]
    duplication_cost=[dup,accumulate_dup+dup]
    #split
    split = split_distance(initial_1,Genome1_block,initial)
    accumulate_split = myTuple[6][1]
    split_cost = [split,accumulate_split+split]
    
    return (initial, elementCount, count,duplication,deletion_cost,duplication_cost,split_cost)
    
''' @function   : find the initial set of blocks of genes/ genes, and provide dictionary that
                  has key is the gene , and value is either 0, 1, 2. 1 means it appears in 1 of them
                  , 2 means it appear in the set and the genome. 
    @input      : 2 tuples of (set,dictionary,Genome count).
    @output     : tuple of (set, dictionary,integer,dictionary))
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

    # create new elementCount dic:
    elementCount = update_dictionary_GG(elementCount1,elementCount2,count)
    ### create the initial Set for the new Tuple:
                
    # edit the initialSet from the Set (reducecount)
    initial_1 = reductionCount(initial1,elementCount)
    #print initial1
    initial_2 = reductionCount(initial2,elementCount)
    # union both of them
    initial = initial_1.union(initial_2)

    ### deal with duplication
    # pull duplication dic from the tuple
    duplication1 = myTuple1[3]
    duplication2 = myTuple2[3]
    # create new duplication dic:
    duplication = update_dictionary_GG(duplication1,duplication2,count)
    if len(duplication)!=0: # only do the reduction Dup if exist duplication 
        initial = reductionDup(initial,duplication)
    initial = reductionSubset(initial,elementCount)

    
    ### calculate edit distances
    #deletion
    deletion = del_distance(initial1,initial2,initial)

    accumulate1_del = myTuple1[4][1] # pull the accumulation del cost
    accumulate2_del = myTuple2[4][1]
    deletion_cost =[deletion,accumulate1_del+accumulate2_del+deletion]
    
    #dup
    if len(duplication)!=0:
        dup = dup_distance(initial1,initial2,initial)
    else:
        dup =0
    accumulate1_dup = myTuple1[5][1]
    accumulate2_dup = myTuple2[5][1]
    duplication_cost=[dup,accumulate1_dup+accumulate2_dup+dup]
    
    #split
    split = split_distance(initial_1,initial_2,initial)    
    accumulate1_split = myTuple1[6][1]
    accumulate2_split = myTuple2[6][1]
    split_cost = [split,split+accumulate1_split+accumulate2_split]
    return (initial, elementCount, count,duplication,deletion_cost,duplication_cost,split_cost)
