# in this version of the class i am removing everyting that has HGT.  This class assumes that the DB searched is protein.

# TODO BLAST field 0 and 1 should actually be generated and formatted by seperate classes that the Homolog class inherits.  In this way any update to information we
# scrape into these fields will automatically be incorporated.  this is not done, and will be put off while i move on, but needs ot be done eventually to make this easier to 
# maintain in the future.


# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment


class Homolog:
    """This is a class that will hold the values that i wish to store about homologs"""
    #def __init__(self, Accession, Organism, Locus, Gene, Predicted_gene, Synonyms, Eval, Percent_ident, Bits_score, GC, Start, Stop, Strand, Product_type, HGT_candidate = {'likelyhood':'not_eval', 'method': 
    #def __init__(self, Accession, Organism, Locus, Gene, Predicted_gene, Synonyms, Eval, Percent_ident, Bits_score, GC, Start, Stop, Strand, Product_type, Alignment_length, Method, Query_accession, Query_common, Query_locus, Query_start): #, Seq):

    def __init__(self, Query_accession, Query_common, Query_locus, BLAST_annotation, Query_start, Query_stop, Query_strand, Query_type, Synonyms, Query_gc, Accession, Organism, Locus, Genbank_annotation, Start, Stop, Strand, GC, Percent_ident, Aligned_length, Number_mismatched, Number_gaps, Align_query_start, Align_query_stop, Align_subject_start, Align_subject_stop, Eval, Bits_score):
        "This will initialize Homolog, the assumed format is all strings, but i am making an  exception for HGT_candidate here since it is an ad-hoc improvement right now"
        
        # information from BLAST field 0: Query
        self.__query_accession = str(Query_accession)
        self.__query_common = str(Query_common)
        self.__query_locus = str(Query_locus)
        self.__blast_annotation = str(BLAST_annotation)
        self.__query_start = int(Query_start)
        self.__query_stop = int(Query_stop)
        self.__query_strand = int(Query_strand)
        self.__query_type = str(Query_type)
        self.__synonyms = Synonyms
        self.__query_gc = float(Query_gc)

        # information from BLAST field 1: subject
        self.__accession = str(Accession)
        self.__organism = str(Organism)
        self.__locus = str(Locus)
        self.__genbank_annotation = str(Genbank_annotation)
        self.__start = int(Start)
        self.__stop = int(Stop)
        self.__strand = int(Strand)
        self.__gc = float(GC)
        
        # The following fields are from the -m 8 fields in order starting at position 2
        self.__percent_ident = float(Percent_ident)
        self.__aligned_length = int(Aligned_length)
        self.__number_mismatched = int(Number_mismatched)
        self.__number_gaps = int(Number_gaps)
        
        self.__align_query_start = int(Align_query_start)
        self.__align_query_stop = int(Align_query_stop)
        self.__align_subject_start = int(Align_subject_start)
        self.__align_subject_stop = int(Align_subject_stop)
        
        self.__eval = float(Eval)
        self.__bits_score = float(Bits_score)
        # currently there are 29 fields

    @classmethod
    def from_file(cls, line):
        try:
            a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,bb = line.strip().split('\t')
            return Homolog(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,bb)
        except:
            print "Error in classmethod from_file in the Homolog class. Number of arguments does not match."

    @classmethod
    def from_blast(cls, line):
        ''''
        print "line.strip().split('\t')", line.strip().split('\t'), len(line.strip().split('\t'))
        query_line, subject_line, percent_ident, aligned_length, number_mismatched, number_gaps, align_query_start, align_query_stop, align_subject_start, \
        align_subject_stop, Eval, bits_score = line.strip().split('\t')
        print "align_query_stop", align_query_stop
        a,b,c,d,e,f,g,h,i,j = query_line.split('|')
        print "subject_line", subject_line, len(subject_line.split('|'))
        print a,b,c,d,e,f,g,h,i,j
        k,l,m,n,o,p,q,r = subject_line.split('|')
        print "k,l,m,n,o,p,q,r", k,l,m,n,o,p,q,r
        return Homolog(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,percent_ident, aligned_length, number_mismatched, number_gaps, align_query_start, \
        align_query_stop, align_subject_start, align_subject_stop, Eval, bits_score)
        '''
        try:
            query_line, subject_line, percent_ident, aligned_length, number_mismatched, number_gaps, align_query_start, align_query_stop, align_subject_start, \
            align_subject_stop, Eval, bits_score = line.strip().split('\t')
            #print "query_line.split('|')", query_line.split('|')
            a,b,c,d,e,f,g,h,i,j = query_line.split('|')
            k,l,m,n,o,p,q,r = subject_line.split('|')
            return Homolog(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,percent_ident, aligned_length, number_mismatched, number_gaps, align_query_start, align_query_stop, align_subject_start, align_subject_stop, Eval, bits_score)
        except:
            print "Error in classmethod from_blast in the Homolog class."
            print line



    # This is information from the query field that is stored in the homolog class
    def query_accession(self):
        return str(self.__query_accession)
        
    def query_common(self):
        return str(self.__query_common)
        
    def query_locus(self):
        return str(self.__query_locus)
    
    def blast_annotation(self):
        return str(self.__blast_annotation)
        
    def query_start(self):
        return int(self.__query_start)
        
    def query_stop(self):
        return int(self.__query_stop)
        
    def query_strand(self):
        return int(self.__query_strand)
        
    def query_type(self):
        return str(self.__query_type)
    
    def synonyms(self): # i should just standardize this outright...
        return self.__synonyms
        
    def query_gc(self):
        return float(self.__query_gc)
    
    
    # This is information from the subject field that is stored in the homolog class                         
    def accession(self):
        """Return the accession number of the organism where the homolog is detected. The return type is string."""
        return str(self.__accession)
    
    def organism(self):
        """Return the organism's common name of the organism where the homolog is detected. The return type is string."""
        return str(self.__organism)
    
    # quick note here, locus must be unique per the definition from NCBI, see the link:
    # http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html    
    def locus(self):
        """Return the locus of the homolog detected. The return type is string."""
        return str(self.__locus)
        
    def genbank_annotation(self):
        return str(self.__genbank_annotation)
        
    def start(self):
        return int(self.__start)
            
    def stop(self):
        return int(self.__stop)
        
    def strand(self):
        return int(self.__strand)
    
    def gc(self):
        return float(self.__gc)
        
    # This returns information about the other BLAST fields relating to the specifics of the hit
    def percent_ident(self):
        return float(self.__percent_ident)
        
    def aligned_length(self):
        return int(self.__aligned_length)
        
    def number_mismatched(self):
        return int(self.__number_mismatched)
        
    def number_gaps(self):
        return int(self.__number_gaps)
    
    def align_query_start(self):
        return int(self.__align_query_start)
    
    def align_query_stop(self):
        return int(self.__align_query_stop)
    
    def align_subject_start(self):
        return int(self.__align_subject_start)
    
    def align_subject_stop(self):
        return int(self.__align_subject_stop)
                    
    def e_val(self):
        return float(self.__eval)

    def bits_score(self):
        return float(self.__bits_score)



    
    def ret_str(self, delim = '\t'):
        return delim.join([str(i) for i in[self.query_accession(), self.query_common(), self.query_locus(), self.blast_annotation(), self.query_start(), self.query_stop(), self.query_strand(), \
            self.query_type(),  self.synonyms(), self.query_gc(), self.accession(), self.organism(), self.locus(), self.genbank_annotation(), self.start(), self.stop(), self.strand(), \
            self.gc(), self.percent_ident(), self.aligned_length(), self.align_query_start(), self.align_query_stop(), self.align_subject_start(), self.align_subject_stop(), \
            self.e_val(), self.bits_score()]])

    # TO DO: overload the print function, so it is not Print.  see http://stackoverflow.com/questions/550470/overload-print-python for details here. 
    def Print(self):
        print self.ret_str()
    
    ######################################################
    # Marked for removal. I don't this it has a purpose. #
    # Test first though. :)                              #
    ###################################################### 
    def ReturnVals(self): # I DO NOT LIKE THIS... 
        return self.query_accession(), self.query_common(), self.query_locus(), self.blast_annotation(), self.query_start(), self.query_stop(), self.query_strand(), self.query_type(),  self.synonyms(), self.query_gc(), self.accession(), self.organism(), self.locus(), self.genbank_annotation(), self.start(), self.stop(), self.strand(), self.gc(), self.percent_ident(), self.aligned_length(), self.align_query_start(), self.align_query_stop(), self.align_subject_start(), self.align_subject_stop(),self.e_val(), self.bits_score()

    # i have no idea what the hell the point of this even is.... i have something above that is capable of dealing with this functionality already.
    def ReturnHomologStr(self, delim = '\t'):
        result = delim.join([self.ReturnVals()])
        return result
        
        
    # This is a function that is to return the contents of the homolog to a file that looks exactly like the BLAST tabular output that the object originated from.
    # It has been tested and validated as of 11/10/2014... what a pain
    def to_file(self):
        
        query_str = '|'.join([str(i) for i in [self.query_accession(), self.query_common(), self.query_locus(), self.blast_annotation(), self.query_start(), self.query_stop(), self.query_strand(), self.query_type(),  self.synonyms(), self.query_gc()]])
        
        subject_str = '|'.join([str(i) for i in [self.accession(), self.organism(), self.locus(), self.genbank_annotation(), self.start(), self.stop(), self.strand(), self.gc()]])
        
        return '\t'.join([str(i) for i in [query_str, subject_str, self.percent_ident(), self.aligned_length(), self.number_mismatched(), self.number_gaps(), self.align_query_start(), self.align_query_stop(), self.align_subject_start(), self.align_subject_stop(), self.e_val(), self.bits_score()]])
        
        
        
        
        
        
        
        
    
        
    
