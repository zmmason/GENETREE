# !/usr/bin/env python3
# University of California, Santa Cruz
# Biomolecular Engineering and Bioinformatics
# Names: Zachary Mason (zmmason), Justin Chan (jumchan)


from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from sequenceAnalysis import OrfFinder, FastAreader
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO, AlignIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# PRELIMINARY PROGRAM

"""
The purpose of the class GeneTree is to have methods to be able to change a list of orfs of
coronavirus sequences to be able to create gene trees shared by the viruses.

The class is made up of three methods, geneSepcificRecord(), fastatoPhylip() and printGeneTree(), and two class
attributes, self.fname and self.newPhylip that is added and moodified to by the three methods to be able to create names
for outout files to be used to create gene trees. The class is based on our previous work of finding open reading frames
in fasta files but, we utlitzed file conversion and phylo in Biopython to create gene trees.

input: A fasta file that contains headers and sequences for coronavirus strains
output: A gene tree or a visual representation of evolutionary relationships of coronaviruses through a gene.

Assumptions:
The coronavirus sequences are DNA sequences or in written in DNA form containing only nucleotides A, T, G and C.
"""
class GeneTree:
    """
    The class, GeneTree, is formatted to have one constructor and three
    methods that format an output file that is able to be aligned to create
    gene tree sequences of corona viruses.
    
    The class contains two attributes, self.fname and self.newPhylip, which are used to be able to read
    and write to them from the three methods in the class. They are names given to files which are assessible
    by each method to either read or modify for its own particular function.
    
    All methods in the classes all have their own purpose where geneSpecificRecord() creates a SeqRecord that fills in gap
    characters for the coronaviruses to be aligned for every gene. fastatoPhylip() uses the SeqRecord to write .phy file that
    will be read by printGeneTree(). printGeneTree() aligns these gene sequences and then constructs evolution distances
    and uses a draw method to construct the gene tree.
    """
    def __init__(self, fname=''):
        """
        The constructor takes in an optional string parameter that represents the name of the .phy
        file.
        
        The constuctor simply is used to set up class attributes and .phy file names so that each method
        can fuifill its own purpose by being able to access these variables.
        """
        self.fname = fname
        nfname = self.fname.split('.', 1)
        newPhylip = nfname[0]+'.phy'  # setting up new file name to be created in directory
        self.newPhylip = newPhylip
        
    def geneSpecificRecord (self, orfList, headList, num):
        """
        Creates a SeqRecord containing aligned sequences with gaps
        
        Using the ORFfinder module from SequenceAnalysis, a target gene can be chosen and identified within the various
        genomes. Data from each gene sequence is compiled into the SeqRecords database after adding gap sequences.
        
        input: An orf list and header list of coronavirus sequences and a number representing which gene is being analyzed
        output: A SeqRecord of a particular gene in coronaviruses containing gap characters 
        
        Keyword Arguments:
        orfList -- A list of a list of orfs in each coronavirus sequence
        headList -- A list of headers of each coronavirus sequence
        num -- The index of a gene in a ordered by descending size of coronaviruses orfs
        
        Assumptions:
        num must be an index of gene that must be present all coronaviruses. In this sample, only the largest four are in
        every sequence.
        """
        sequenceInfo = []
        for gene in orfList:  # Finds target gene in each genome
            sequenceInfo.append(gene[num]) # ***any gene can be utilized***
        longestLength = max(len(s) for s in sequenceInfo) # gets longest seq to match length with gap characters
        paddedSequences = [s.ljust(longestLength, '-') for s in sequenceInfo]  # Adds gap characters
        
        records = (SeqRecord(Seq(s), id = str(paddedSequences.index(s))) for s in paddedSequences)  #creating a SeqRecord
        return(records)
    
    def fastaToPhylip(self, records):
        """
        Converts and aligns sequences from .fasta formatting to .phy formatting
        
        Writes to an output .phy file using a list of Seq objects called records. Similar to create a stdout file later
        used to construct gene trees in printGeneTree()
        
        input: A SeqRecord of coronaviruses
        output: A .phy file containing a SeqRecord
        
        Keyword Arguments:
        records -- SeqRecord contain Seq objects that each represent a coronavirus
        
        Assumptions:
        None
        """
        SeqIO.write(records, self.newPhylip, 'phylip') # Writes a new .phy file containing the SeqRecord for specific gene
    
    def printGeneTree(self):
        """
        Prints out gene trees with matplotlib and in the terminal for the four largest ORFs of coronaviruses
        
        Takes a .phy file containing multiple alligned sequences, generates a matrix based on sequence composition 
        and compares each sequence (genome) to one another. sequences with grater scores (similarity) are ranked closer
        together on the phylogenetic trees.
        
        input: A .phy file that contains coronavirus gene sequences to draw phylogenetic tree
        output: A visual representation of a gene tree on terminal and matplotlib
        
        Assumptions:
        The .phy file was formatted correctly with the right genes
        """
        align = AlignIO.read(self.newPhylip, 'phylip')  # Reads created .phy file containing the SeqRecord
        #print (align) # prints concatenated allignments
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(align)# Calculate the distance matrix
        print('\n======================================== DISTANCE MATRIX =======================================\n')
        print(dm,"\n\n") # Print the distance Matrix
        constructor = DistanceTreeConstructor() # Construct the phylogenetic tree using UPGMA algorithm
        tree = constructor.upgma(dm)
        print('\n========================================= GENE TREE ===========================================\n')
        Phylo.draw(tree) # Draw the phylogenetic tree (must install matplotlib to use this formatting)
        Phylo.draw_ascii(tree) # Print the phylogenetic tree in terminal
        
###############################################################################################################################
# Main
#
# The purpose of the main method is to format the output and call methods to be able to set up the process of drawing the gene
# trees for the four largest open reading frames that all corona virus in the file contain by creating a FastAreader() and
# GeneTree().
#
# FastAreader is used to parse the fasta file which contains the seven coronaviruses combined into headers and sequences
# which are used for creating the key and gene trees. The sequence is then analyzed for open reading frames and then used
# to create SeqRecord and then phy. file that is then aligned through a multiple sequence alignment to be able to draw 
# the gene trees for the four largest open reading frames: ORF1A, ORF1B, ORF S, ORF N
#
############################################################################################################################### 

def main(inCL=None):
    """
    Uses a fasta file to create gene trees for the four largest genes on coronaviruses
    
    This program takes a "combined" fasta file containing multiple genomes from similar virus strains and searches each 
    genome for a common specific gene. These genes are then aligned to one another and converted to a .phy file that can
    then be used to create a detailed phylogenetic tree based on the variation between particular genes. In theory, the
    variation in the molecular composition of the genes will determine the trend in the variation of the genomes. This allows
    the "Gene Tree" to provide a rough outline of what the "Phylogenic Tree" would look like having alligned full genomes.
    
    input: An input fasta file that contains virus headers and their DNA sequence
    output: A written and graphical description of four gene trees of the four largest orfs in coronaviruses. The output also
    contains a distance matrix for each gene tree and key correlating numbers to viruses.
    
    Assumptions:
    - input file must follow fasta format and must be a DNA nucleotide sequence.
    - input file contains solely corona viruses since they contain only 4 genes that are
    conserved among each other
    """

    headList = []  # Stores header of coronavirus sequences in fasta file
    orfList = []  # Stores sequences containing ORFs of coronavirus sequences in fasta file
    validNucs = ['A', 'C', 'G', 'T']
    myReader = FastAreader('Combined-ALL-SARS-CoV.fasta') 
    for head, seq in myReader.readFasta(): # Using fastAreader to read in .fasta files
        headList.append(head)
        for i in seq:
            if i not in validNucs:  # Removing non-valid bases
                seq = seq.replace(i,"") 
        orf = OrfFinder(seq, 300, True)  # Includes the largest ORF greater than 300 nucleotides within a stop codon
        geneOrfList = orf.getOrfs()
        geneSeq = []  # Stores ORF sequences
        for openFrame in geneOrfList:
            geneSeq.append(seq[openFrame[1]-1:openFrame[2]-1])
        orfList.append(geneSeq)
     
    # Calls methods to create SeqRecords and then .py file to print gene trees
    myPhylo = GeneTree() 
    for i in range(0,4,1):  # Loops to print the first four gene trees of every sequence
        records = myPhylo.geneSpecificRecord(orfList, headList, i)  # Creates list of SeqRecords that represent a sequence
        # alignments = myPhylo.fastaToPhylip(records)  # Makes a .phy file using a .fasta file
        print("GENE " + str(i+1) + ":")
        # printTree = myPhylo.printGeneTree()  # Prints Gene Trees
    
    x = 0
    print('\n\n============================================ K E Y ============================================\n')
    for header in headList: # Loops through headers to print key
        header = header.split(',')
        header = header[0]
        print("{} = {}" .format(x, header))  # Prints each line containing the header
        x += 1


if __name__ == "__main__":
    main()  # Runs the main method to print output with gene trees