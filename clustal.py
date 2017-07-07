############################################
# Clustal Omega Alignment Parser           #
# Written by Hasan Alkhairo & Skyler Kuhn  #
# VCU Bioinformatics MS Candidates         #
# Version 1.1.0 Released 01/19/2017 Spring #
############################################
### Imports ###
from __future__ import division, print_function
import re
import subprocess as sp
import sys

### Objects & Methods ###
class ClustalReport(object):
    """
    This creates a Clustal Report object which takes in a clustal alignment file
    then parses it extracts relevant information
    """
    def __init__(self,alnfile):
        self.alnfile = alnfile
        self.alignments = []
        self.headers = []
        self.score_raw = ''
        self.parse(alnfile)
        self.type = self.find_type()

    def __len__(self):
        return len(self.alignments)

    def find_type(self):
        if len(self.alignments) > 0:
            if re.search('^[ACGTN\-]+$',self.alignments[0]):
                return "DNA"
            elif re.search('^[ARNDCEQGHILKMFPSTWYVOUXB\-]+$',self.alignments[0]): #OUXB are non standard AAs
                return "Amino Acids"
        return "Unknown"

    def score(self,ignore_gaps = False):
        if self.type == "DNA":
            score_dict = {'*' : 10,':' : 0, '.' : 0 }
        else:
            score_dict = {'*' : 10,':' : 5, '.' : 1 }
        score_len = 0
        score = 0
        for char in self.score_raw:
            score += score_dict.get(char,0) #sums all scores from each position
            if ignore_gaps:
                if re.search('\s',char):
                    continue
                else:
                    score_len += 1
            else:
                score_len += 1
        return score/score_len #returns average

    def parse(self,alignment = None):
        if not alignment:
            alignment = self.alnfile
        for line in open(alignment):
            if re.search("CLUSTAL O",line,re.I):
                continue
            else:
                line = line.strip()
                if re.search('^\s*[\.\:\*][^A-Za-z]\s*',line):
                                        self.score_raw = line #score_raw will be reassigned here to the last line of each aln file, which includes the raw scores.
                                        #The raw scores will be converted to assigned numbers and an average will be calculated in the Score function.
                else:
                    line = re.split('\s+',line)
                    if len(line) > 1:
                        self.headers.append(line[0])
                        self.alignments.append(line[1])

class ClustalPipeline(object):
    """
    This is designed to validate an input fasta file (not implemented yet)
    Then execute clustal omega to align the sequences within the file
    """
    def __init__(self, infile, outfile = 'out.aln', threads = 4):
        self.infile = infile
        self.outfile = outfile
        self.threads = threads
    def validate(self): #come back to and delete
        pass

    def run(self):
        self.validate()
        #callling subprocess to run clustal
        #Usage: list of arguments
        sp.Popen(['/usr/bin/clustalo','-i',self.infile,'--threads=%i' % self.threads,'-o',self.outfile,'--outfmt=clustal','-v','--force','--wrap=999999']).wait()

### Main Method if script is being used independently ###
if __name__ == '__main__':
        try:
                infile = sys.argv[1]
                outfile = sys.argv[2]
        except:
                print ("No infile and outfile given, invoking demo mode... using bat_flu.fasta & testout.aln")
                print ("Next time specify a file to be aligned & and output filename")
                infile = 'Segment_8_seq.fasta'
                outfile = 'testout.aln'
        print("Intializing and running Clustal Omega Pipeline")
        clustalpipe = ClustalPipeline(infile,outfile,10)
        clustalpipe.run()
        print("Finished Running Clustal Pipeline.")
        print("Beginning Parsing")
        alignment_report = ClustalReport(outfile)
        print ("Found alignment of: ",alignment_report.type,"With: ",len(alignment_report),"alignments.")
        print ("Overall alignment score: ",alignment_report.score())
        print ("Done")