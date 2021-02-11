#!/usr/bin/env python

import argparse
import re
import cairo
import math
import seaborn as sns


def get_args():
    '''to run script from the command line'''
    parser = argparse.ArgumentParser(description="find motifs and exons in fasta file, outputs a diagram of motif/exon locations")
    parser.add_argument("-f", "--fasta", help="input fasta file")
    parser.add_argument("-m", "--motifs", help="txt file of motifs")
    return parser.parse_args()

args = get_args()

input_file = args.fasta
motif_file = args.motifs

output_file = input_file.split(".")[0] + ".svg"


def DegNucleoConverter(motif):
    '''Take a motif as input and returns new motif containing regex characters for any degenerate nucleotides'''
    IUPAC_dict = {'U': '[TU]', 'T': '[TU]', 'W': '[ATU]', 'S': '[CG]', 'M': '[AC]', 'K': '[GTU]', 'R': '[AG]', 'Y': '[CTU]', 'B': '[CGTU]', 'D': '[AGTU]', 'H': '[ACTU]', 'V': '[ACG]', 'N': '[ACGTU]'}
    motif = motif.upper()
    motif_split = list(motif)
    for n, i in enumerate(motif_split):
        if i in IUPAC_dict:
            motif_split[n] = IUPAC_dict[i]
    motif = ''.join(motif_split)
    return motif

def MotifFileParser(motif_file):
    '''Takes a motif file as input and returns a dictionary of converted motifs with their original degenerate sequences as values'''
    with open (motif_file, "r") as mf:
        motif_dict = {}
        num_motifs = 0
        for line in mf:
            motif = line.strip()
            num_motifs += 1
            conv_motif = DegNucleoConverter(motif).lower()
            motif_dict.setdefault(conv_motif, motif)
    return motif_dict, num_motifs

motif_dict, num_motifs = MotifFileParser(motif_file)

def MotifFinder(sequence):
    '''Find start locations of motifs in a read and returns dictionary of start indexes'''
    motif_locations = {}
    for motif in motif_dict:
        locations = [i.start() for i in re.finditer(motif, sequence)]
        motif_locations.setdefault(motif_dict[motif], locations)
    return motif_locations

def ExonFinder(header, sequence):
    '''Takes a sequences as input and returns index of exon and exon length'''
    gene_name = (header.split(" ")[0])[1:]
    exon_location = []
    read_length = 0
    for n, i in enumerate(sequence):
        read_length += 1
        if i.isupper():
            exon_location.append(n)
    #creates tuple containing location of start of exon and exon length
    exon_info = (exon_location[0], len(exon_location), gene_name)
    return exon_info, read_length

def DrawGene(read_length, exon_info, motifs, read_counter, end_of_file, surface, context):
    '''Plots gene with motifs and exon using pycairo'''
    #grab gene name
    gene_name = exon_info[2]
    #create color palette for motifs and title legend
    motif_color_pal = sns.color_palette("hls", num_motifs)
    context.select_font_face("Sans", cairo.FONT_SLANT_NORMAL,
                        cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(16)
    context.move_to(50, 20)
    context.set_source_rgba(0, 0, 0)
    context.show_text("Motif key:")
    #write gene name
    x, y = 50, 150+(80*read_counter)

    context.set_font_size(14)
    context.move_to(x, y-20)
    context.show_text(gene_name)
    #draw gene as line and initialize starting point
    context.set_line_width(5)
    context.move_to(x, y)
    context.line_to(x+read_length, y)
    context.stroke()
    #draw exon as thick line
    context.move_to(x+exon_info[0], y)
    context.set_line_width(30)
    context.line_to(x+exon_info[0]+exon_info[1], y)
    context.stroke()
    #draw motifs
    motif_counter = 1
    for i in motifs:
        x1, y1, y2 = x, y+15, y-15
        context.set_line_width(len(i))
        motif_counter += 1
        motif_color = motif_color_pal[motif_counter-2]
        context.set_source_rgba(motif_color[0], motif_color[1], motif_color[2])
        for j in motifs[i]:
            context.move_to(x+j, y1)
            context.line_to(x+j, y2)
            context.stroke()
        #create motif key at top of page
        x0, y0 = 50, 25*motif_counter
        context.set_line_width(20)
        context.move_to(x0, y0)
        context.set_font_size(12)
        context.line_to(x0+20, y0)
        context.stroke()
        context.move_to(x0 + 30, y0)
        context.set_source_rgba(0,0,0)
        context.show_text(i.upper())
    #write out surface if end of file
    if end_of_file == True:
        surface.finish()
    return surface

def FastaParser(input_file, output_file):
    '''Create svg output file from input fasta file'''
    records = {}
    read_counter = 0
    surface = cairo.SVGSurface(output_file, 1000, 1000)
    #create starting surface
    context = cairo.Context(surface)
    with open(input_file) as fh:
        for line in fh:
            end_of_file = False
            line = line.strip()
            #check if line is a header line
            if line[0] == '>':
                #check to see if header variable exists yet
                if "header" in locals():
                    read_counter += 1
                    #find motifs in current record
                    motif_locations = MotifFinder(records[header].lower())
                    #find exon location
                    exon_info, read_length = ExonFinder(header, records[header])
                    #draw gene with exons and motifs
                    DrawGene(read_length, exon_info, motif_locations, read_counter, end_of_file, surface, context)
                header = line
                records.setdefault(header, "")
            else:
                records[header] += str(line)
        #used at end of file for last record        
        else:
            read_counter += 1
            end_of_file = True
            #find motifs in current record
            motif_locations = MotifFinder(records[header].lower())
            #find exon locations
            exon_info, read_length = ExonFinder(header, records[header])
            #draw gene with exons and motifs
            DrawGene(read_length, exon_info, motif_locations, read_counter, end_of_file, surface, context)
    return motif_locations

motif_locations = FastaParser(input_file, output_file)

            
            
