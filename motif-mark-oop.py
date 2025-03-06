#!/usr/bin/env python

import argparse     #for command line options
import cairo        #to draw output image
import re           #regex


##########################
#        ARGUMENTS       #
##########################


def get_arg():
   parser = argparse.ArgumentParser(description="")
   parser.add_argument("-f", help="designates absolute file path to fasta file", required=True, type=str)
   parser.add_argument("-m", help="designates absolute file path to motif file", required=True, type=str)
   return parser.parse_args()

args = get_arg()
f=args.f
m=args.m


##########################
#    GLOBAL VARIABLES    #
##########################


FEATURE_HEIGHT = 25     #height of exons on image
CENTERED = 50       #draws a gene line every 50

COLORS = [(0/255, 192/255, 0/255), (80/255, 208/255, 255/255), (160/255, 32/255, 255/255), (255/255, 96/255, 208/255), (255/255, 160/255, 16/255)]      #colors for up to 5 unique motifs, script must be capable of handling max 5 motifs

IUPAC = {
    "U":"[UT]",
    "W":"[ATU]",
    "S":"[CG]",
    "M":"[AC]",
    "K":"[GTU]",
    "R":"[AG]",
    "Y":"[CTU]",
    "B":"[CGTU]",
    "D":"[AGTU]",
    "H":"[ACTU]",
    "V":"[ACG]",
    "N":"[ACGTU]",
}
#creating lists of all possible nucleotides associated with IUPAC ambiguity characters

##########################
#         CLASSES        #
##########################


class Gene:

    def __init__(self, start, end, gene_name, gene_num):        #class Gene requires a start position, end position, gene name (for labeling), and a gene number to iterate through the fasta file
        '''Motif sequences'''

        self.start = start
        self.end = end
        self.gene_name = gene_name
        self.gene_num = gene_num

#drawing lines of genes, proportional in size to the fasta length
    def draw(self, context):
        context.set_line_width(1)
        context.set_source_rgb(0, 0, 0)
        y = (self.gene_num - 1)*100 + CENTERED*2
        context.move_to((self.start + CENTERED), y)
        context.line_to(self.end+CENTERED,y)
        context.stroke()

#labeling gene lines with gene names
        context.select_font_face("Verdana")
        context.set_font_size(10)
        context.move_to(25, y - 30)
        context.show_text(self.gene_name)

#adding title to figure, doesn't have to be in this class but it is
        context.select_font_face("Verdana")
        context.set_font_size(15)
        context.move_to(25, 15)
        context.show_text("Motif Mark Visualization of " + args.f) #title of image will depend on input file


class Exon:

    def __init__(self, start, end, gene_num):       #class Exon requires a exon start, exon end, and gene number for iterating
        '''Motif sequences'''

        self.start = start
        self.end = end
        self.gene_num = gene_num

#drawing exons on genes
    def draw(self, context):
        context.set_source_rgb(0, 0, 0)
        y = ((self.gene_num - 1) * 100) + CENTERED/2
        context.rectangle(self.start + CENTERED, y+CENTERED*1.25, self.end - self.start, FEATURE_HEIGHT)        #1.25 was found through trial and error
        context.fill()


class Motif:

    def __init__(self, start, end, gene_num, color):        #class Motif requires start of motif, end of motif, gene number for iterating, and color to distiguish between motifs
        '''Motif sequences'''

        self.start = start
        self.end = end
        self.gene_num = gene_num
        self.color = color

#drawing motifs on genes
    def draw(self, context):
        context.set_source_rgb(self.color[0], self.color[1], self.color[2])
        y = ((self.gene_num - 1) * 100) + CENTERED/2
        context.rectangle(self.start + CENTERED, y+CENTERED*1.25, self.end - self.start, FEATURE_HEIGHT)
        context.fill()


##########################
#        FUNCTIONS       #
##########################


def find_exon(sequence: str) -> tuple:
    '''Finds exons in fasta sequences to map on genes. Exons are distiguished in fasta sequences with capitalization.'''

#for each index in the length of the fasta record, if the character is upper case save that index as 'start'
    for i in range(len(sequence)):
        if sequence[i].isupper():
            start = int(i)
            break

#as soon as the character becomes lower case, save the index as 'end' of the exon and break the loop
    for i in range(start, len(sequence)):
        if sequence[i].islower():
            end = int(i)
            break

    return (start,end)


def possible_motifs(motif:str) -> str:
    '''Finds motifs in fasta sequences and replaces ambiguous nucleotides with all possible IUPAC options.'''

    motif = motif.upper()
    for ambiguous_nuc, iupac_options in IUPAC.items():      #for each ambiguous base and its corresponding iupac options in the iupac dictionary
        motif = motif.replace(ambiguous_nuc, iupac_options)     #replace the ambiguous base with all iupac options for that base
    
    return(motif)


##########################
#          MAIN          #
##########################


fasta_dict:dict = {}
with open (f, "r") as fr:       #opening fasta file for reading
    while True:
        line = fr.readline().strip()        #for each line in the fasta file, strip the newline character at the end

        if line == "":      #break when reach end of file
            break

        #saving header line of each fasta record as the key in the dictionary
        if line.startswith(">"):
            header = line
            fasta_dict[header] = ""

        #saving the sequence as the value for the key in the dictionary
        else:
            fasta_dict[header] += line


motif_dict = {}
with open(m, "r") as mf:
    while True:
        line = mf.readline().strip()

        if line == "":
            break
        
        #motif file should only hold one motif per line
        #real motif is the key, iupac motif options are the value
        motif_dict[line] = possible_motifs(line)


#finding length of longest seqeunce to determine width of cairo canvas
long_seq = 0

#iterate through each record until the longest sequence length is set as 'long_seq'
for key in fasta_dict:
    if len(fasta_dict[key]) > long_seq:
        long_seq = len(fasta_dict[key])

#initializing cairo figure
width = long_seq + 225      #width according to length of longest sequence plus buffer room
height = 100 * len(fasta_dict) + 50     #adding 100 units for each fasta record present plus buffer room
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32,width, height)
context = cairo.Context(surface)
context.set_source_rgb(1, 1, 1)     #making background white
context.paint()

#finding genes in fasta_dict and drawing according to class Gene
for i, header in enumerate(fasta_dict):
        gene = Gene(0, len(fasta_dict[header]), header, i+1)
        gene.draw(context)

        #find start and end of exon, drawing according to class Exon
        (start, end) = find_exon(fasta_dict[header])
        exon = Exon(start, end, i+1)
        exon.draw(context)

        #find motifs in sequences, draw according to class Motif
        for j, key in enumerate(motif_dict):
            matches = re.finditer(motif_dict[key], fasta_dict[header].upper())
            for match in matches:
                start = match.start()
                end = match.start() + len(match.group())
                motif = Motif(start, end, i+1, COLORS[j])
                motif.draw(context)


##########################
#           KEY          #
##########################

#title for key
context.select_font_face("Verdana")
context.set_font_size(15)
context.set_source_rgb(0, 0, 0)
context.move_to(width-125, 15)
context.show_text("Key")


#drawing out colored boxes and text for motif key
for i, motif in enumerate(motif_dict):
    context.set_source_rgb(COLORS[i][0], COLORS[i][1], COLORS[i][2])
    y = ((i + 1) * 30) + 10
    context.move_to(width-100, y)
    context.show_text(motif)
    context.rectangle(width-125, y-10, 10, 10)
    context.fill()


name = f.split(".")     #naming .png file based on input fasta file
surface.write_to_png(f"{name[0]}.png")





