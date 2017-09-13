#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Automate create group dic for the chosen taxa, write it to a file
              I classift by class
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
    myclass=set() # keep track of class
    class_dic={} # key is accesion number, value is the class
    # color avaliable for SVG_Color of ete3
    color_list=['hotpink','deepskyblue','black','brown','yellow','magenta','purple',
               'green','mediumblue','silver']
    color_dic={}
    if condition:
        accession = parse_accession(args.AccessionNumber)
        outputsession = args.OutputFile
        res = traverseAll(args.InputGenBackDirectory)
        index =0
        for r in res:
            root,f = os.path.split(r)
            accession_num= f.split('.')[0]
            # check if the file name is in our ancession file
            if accession_num in accession: # NC_000964.gbk
                ''' parse using SeqIO on genbank is too slow. I will open
                    the file and read it directly 
                genome=SeqIO.parse(r,'genbank')
                first_rec = next (genome)
                genome_class = first_rec.annotations['taxonomy'][2]'''
                
                # dirty version of manipulate genbank file
                myfile = open(r,'r')
                flag= False
                for line in myfile.readlines():
                    if flag:
                        info = line
                        break
                    if 'ORGANISM' in line:
                        flag= True
                genome_class =str(info.strip(' ').split(';')[2])
                class_dic[accession_num]=genome_class
                myclass.add(genome_class)
        # assign the collor to the class
        for item in myclass:
            # randomly choose a color
            # get the color
            color_dic[item]=color_list[index]
           # remove the color
            color_list.remove(color_list[index])
            
        # assign the accession the the color:
        for key in class_dic:
            color_dic[key]=color_dic[class_dic[key]]
        # writing to outfile
        outfile=open(outputsession,'w')
        for key in color_dic:
            string = key +':'+ color_dic[key]+'\n'
            outfile.write(string)
        outfile.close()
