#!/usr/bin/env python
import argparse
import itertools
import re
import gzip
import csv
from argparse import RawTextHelpFormatter

def get_args():
    parser = argparse.ArgumentParser(description=str("This program returns a .svg file portraying motif locations along introns and exons for a specific gene"), formatter_class=RawTextHelpFormatter)
    parser.add_argument("-f", "--FASTAfile", help="FASTA file containing genes to be examined.", required=True)
    parser.add_argument("-m", "--motiffile", help="File containing known motifs. No header information, each motif on its own line.", required=True)

    return parser.parse_args()


args = get_args()
input_file = args.FASTAfile
motif = args.motiffile


def parse_fasta():
    with open(input_file, "r") as f:
        counter = 0
        buffer1 = []
        buffer2 = []
        headers = []
        for line in f:

            #Define first use case because first if statement is "not in". If first use not defined, sequences from the first record will be recorded properly
            counter += 1
            if counter == 1:
                line = line.strip()
                headers.append(line) #Append the gene name from the first line to headers
                continue

            if ">" not in line:
                line = line.strip()
                buffer1.append(line) #Collect seqs
            else:
                line = line.strip()
                cleanbuff = ""
                for i in buffer1:
                    cleanbuff += str(i) # Concatenate all lines from same fasta record as a string
                buffer2.append(cleanbuff) #Append all protein seqs from previous record to one list element as a string
                headers.append(line) #Append header to list
                buffer1 = [] #Clear buffer1
        cleanbuff = ""
        for i in buffer1:
            i = i.strip()
            cleanbuff += str(i)
        buffer2.append(cleanbuff)
        
        fasta_dict = {}
        for i in range(len(headers)):
            fasta_dict[headers[i]] = buffer2[i]

    return fasta_dict

def parse_motif():
    '''This function takes a motif file and parses each individual line into its own list item in a list called motifs'''
    with open(motif, "r") as m:
        motifs = []
        for line in m:
            line = line.strip()
            motifs.append(line)
    return motifs

def translate_motifs():

    '''This function takes a list of motifs (one per list item) and translates it according to the IUPAC definitions. 
    This function outputs a list called translated_motif_regex with items that can be put into a regular expression find function. 
    Example: "Y" indicates the base can either be cytosine or thymine. It must be one or the other. It can be capitalized or not. 
    Therefore, "Y" = [cCtTuU] in a regex function'''

    IUPAC_dict = {    
    "A":"[Aa]",
    "C":"[Cc]",
    "G":"[Gg]",
    "T":"[tTuU]",
    "U":"[uUtT]",
    "W":"[aAtTuU]",
    "S":"[cCgG]",
    "M":"[aAcC]",
    "K":"[gGtTuU]",
    "R":"[aAgG]",
    "Y":"[cCtTuU]",
    "B":"[cCgGtTuU]",
    "D":"[aAgGtTuU]",
    "H":"[aAcCtTuU]",
    "V":"[aAcCgG]",
    "N":"[aAcCgGtTuU]",
    "Z":"[-]",
    }

    translated_motif_regex = []
    for item in motifs:
        item = item.upper()
        regex_motif = []
        regex = ""
        for letter in item:
            regex_motif += IUPAC_dict[letter]
        
        for item in regex_motif:
            regex += str(item)
        
        translated_motif_regex.append(regex)

    return translated_motif_regex


def find_exons():
    '''This function finds the start and stop locations of exons within each gene and puts that information as values in a dictionary where gene names are the key.'''
    exon_dict = {}
    for i in fasta_dict:
        x = fasta_dict[i]
        iterator = re.finditer("[A-Z]+", x)
        for match in iterator:
            exon_dict[i] = match.span()
    
    return exon_dict

def find_longest_sequence():
    '''This function finds the longest gene with which to scale the graphic off of for all genes. Each pixel will be 1 nucleotide.'''
    longest_gene = 0
    for i in fasta_dict:
        x = fasta_dict[i]
        if len(x) > longest_gene:
            longest_gene = len(x)
    
    return longest_gene



############# MAIN FUNCTION ##################

fasta_dict = parse_fasta()

motifs = parse_motif()

translated_motif_regex = translate_motifs()

exon_dict = find_exons()

longest_gene = find_longest_sequence()


import cairo
colors_list = [0,0,1,0,1,0,1,0,0,1,0,1,1,1,0,0,1,1] # Set initial colors to every permutation of 1,0,0 and 1,1,0 for best differentiation of colors
counter = 0
if len(translated_motif_regex) > 6: #If there are more than 6 motifs, add more colors to the color list
    for i in range(len(translated_motif_regex) - 6):
        x = random.randrange(3,10,1)/10

        if counter % 3 == 0:
            colors_list.append(x)
            colors_list.append(x)
            colors_list.append(0)
        elif counter % 3 == 1:
            colors_list.append(0)
            colors_list.append(x)
            colors_list.append(x)
        elif counter % 3 == 2:
            colors_list.append(x)
            colors_list.append(0)
            colors_list.append(x)
        counter += 1


color_counter = 0
width = longest_gene + 100
height = len(fasta_dict) * 100 + len(translated_motif_regex) * 25 + 75
fasta_count = 0

short_name = re.search("^(.*)\.", input_file)
short_name = short_name.group(1)

surface = cairo.SVGSurface("./{fname}.svg".format(fname = short_name), width, height)


context = cairo.Context(surface)

current_height = 50

for i in fasta_dict:
    if fasta_count > 0: # restore black color fill
        context.restore()
    a = fasta_dict[i]
    context.move_to(50, current_height + 25) #Set context to draw gene line
    context.set_line_width(5) #Set width of gene line
    context.line_to(len(a) + 50,current_height + 25) # Draw Gene Line
    context.stroke()
    context.move_to(50, current_height + 25 - 40)
    context.show_text(i)
        
    context.set_line_width(1) #Set exon box line width
    exon_coordinates = exon_dict[i]
    list(exon_coordinates)
    context.rectangle(exon_coordinates[0] + 50, current_height + 25 - 25, exon_coordinates[1] - exon_coordinates[0], 50) #Draw exon box
    context.fill()
    context.save() # Save black color fill for genes and exons
        
    fasta_count +=1

    for i in translated_motif_regex:
        iterator = re.finditer(i, a) #Find matches of each translated motif in each sequence
        for match in iterator:
            match_list = match.span() # return coordinates of each match for each motif in each sequence
            list(match_list)
                
            context.set_line_width(.1)
            context.rectangle(match_list[0] + 50, current_height + 25 - 30, match_list[1] - match_list[0], 60) # Draw rectangle of width(motif) at each motif location
            context.set_source_rgba(colors_list[color_counter], colors_list[color_counter + 1], colors_list[color_counter + 2], 0.65) #Draw in different colors for each motif
            context.fill()
        color_counter += 3
    current_height += 100
    color_counter = 0
    
    
#LEGEND CREATION
color_counter = 0
legend_height = height - len(translated_motif_regex)*50 + 75#Start legend figure below the last gene displayed
context.set_source_rgb(0,0,0) #Set legend text and lines to black
    
context.move_to(50, legend_height - 20)
context.show_text("Legend (motif name x color)") #Title the legend
context.set_line_width(3)
context.move_to(50, legend_height - 10)
context.line_to(width, legend_height- 10)
context.stroke()
    
for i in motifs: #Create legend figure
    context.restore() #Restore black fill for text
        
    context.move_to(50, legend_height + 2.5)
    context.show_text(i) #Display motif name
        
    context.move_to(50, legend_height - 10)
    context.line_to(width, legend_height - 10)
    context.stroke()
    context.save()
        
    context.rectangle(width - 50, legend_height - 5, 10, 10) #Display color of motif lines
    context.set_source_rgb(colors_list[color_counter], colors_list[color_counter + 1], colors_list[color_counter + 2])
    context.fill()
        

        
    color_counter +=3
    legend_height = legend_height + 20
        
        
    
    
surface.finish()