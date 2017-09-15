#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Automate create group dic for the chosen taxa, write it to a file
              I classify them by class
    Start   : 05/10/2016
    End     : 05//2016
'''
import os
import argparse
import uuid
from Bio import SeqIO
# traverse and get the file
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+'/'+f)
    return res

class readable_dir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
           try:
               os.mkdir(prospective_dir)
           except OSError:
               print (argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir)))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--InputGenBackDirectory","-i",action=readable_dir,help="Genbank directory")
    parser.add_argument("--OutputFile","-o", help="Output of this program will be stored in the path supplied here. It will make a new directory if path given is valid or it will raise an error")
    parser.add_argument("--AccessionNumber","-a", help="Accession name (phylo_order.txt)")
    args = parser.parse_args()
    return args


def chk_output_directory_path(OutputDirectory,sessionID):
    if not os.path.exists(OutputDirectory + "_" + str(sessionID)):
        try:
           #os.mkdir(OutputDirectory + "_" + str(sessionID))
           return True
        except OSError:
           print ("Unable to create file:", OutputFile)
           sys.exit()
           
def parse_accession(myFile):
    accession = open(myFile,'r')
    filter_acession = []
    for line in accession.readlines():
        accession_number = line.split('\n')[0]
        filter_acession.append(accession_number)
    accession.close()
    return filter_acession
    
if __name__ == "__main__":

    args = get_arguments()
    sessionID = uuid.uuid1()
    condition = chk_output_directory_path(args.OutputFile,sessionID)
    accessions=args.AccessionNumber
    myclass={0:set(),1:set(),2:set(),3:set(),4:set(),5:set(),6:set(),7:set()} # keep track of class
    class_dic={} # key is accesion number, value is the class
    # color avaliable for SVG_Color of ete3
    color_list=['hotpink','deepskyblue','black','brown','yellow','magenta','purple',
               'green','mediumblue','silver']
    color_dic={}
    if condition:
        accession = parse_accession(args.AccessionNumber)
        outputsession = args.OutputFile
        res = traverseAll(args.InputGenBackDirectory)

        for r in res:
            print (r)
            input_seq_iterator = SeqIO.parse(r, "genbank")
            first_rec = input_seq_iterator.next()
            accession_num= first_rec.annotations['accessions'][0]
            # check if the file name is in our ancession file
            if accession_num in accession: # NC_000964.gbk
                print (accession_num)
                taxonomy = first_rec.annotations['taxonomy']
                index = 0
                for item in taxonomy:
                    myclass[index].add(item)
                    index+=1
                class_dic[accession_num]=taxonomy
        # assign the collor to the class
        for group in sorted(myclass):
            if len(myclass[group])>=4:
                break
        print (myclass)
        index =0
        for item in myclass[group]:
            # randomly choose a color
            # get the color
            color_dic[item]=color_list[index]
           # remove the color
            color_list.remove(color_list[index])
        print (color_dic)
        # assign the accession the the color:
        for key in class_dic:
            color_dic[key]=color_dic[class_dic[key][group]]
        # writing to outfile
        outfile=open(outputsession,'w')
        for key in color_dic:
            string = key +':'+ color_dic[key]+'\n'
            outfile.write(string)
        outfile.close()
