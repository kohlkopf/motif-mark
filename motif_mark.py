#!/usr/local/bin python

import re
import cairo
import math
import numpy as np
import argparse
import webcolors
import randomcolor
import textwrap

def get_arguments():

    parser = argparse.ArgumentParser(description="Creates a visualization of one or more motifs across a sequence or multiple sequences provided in UCSC format. Requires sequences and a list of motifs to be visualized.")
    parser.add_argument("-f", help="File path of asta file to be processed. <str>", required=True, type=str)
    parser.add_argument("-m", help="File path of motif file to define motif sequences to search for. <str>", required=True, type=str)

    return parser.parse_args()

class Sequence():
    '''Object with sequence lines per single fasta ID. Lines merged in to single string for searching.'''
    
    def __init__(self, lines):
        self.gene = lines[0].split(' ')[0].strip('>')
        self.seq = ''.join([line.strip('\n') for line in lines if line.startswith('>') != True])
    
    def tot_seq_len(self):
        return int(len(self.seq))
    
    def exon_bounds(self):
        exon_st = self.seq.find(re.search('[ATCG]+',self.seq)[0])
        exon_ed = (self.seq[exon_st:].find(re.search('[atcg]+',self.seq[exon_st:])[0]))+ exon_st
        return [exon_st, exon_ed]
    
    def exon_bounds_rel(self, trimmed_seq):
        exon_st = trimmed_seq.find(re.search('[ATCG]+',trimmed_seq)[0])
        exon_ed = (trimmed_seq[exon_st:].find(re.search('[atcg]+',trimmed_seq[exon_st:])[0]))+ exon_st
        return [exon_st, exon_ed]
        
def iupac_trans(motifs_file):
    '''Returns list of regular expression terms for the motifs in motif_list.txt based on IUPAC.'''
    
    iupac = {"A":'[Aa]', "a":'[Aa]',
             "C":'[Cc]', "c":'[Cc]',
             "G":'[Gg]', "g":'[Gg]',
             "T":'[TUtu]', "t":'[TUtu]',
             "U":'[TUtu]', "u":'[TUtu]',
             "R":'[AGag]', "r":'[AGag]',
             "Y":'[CTct]', "y":'[CTct]',
             "S":'[GCgc]', "s":'[GCgc]',
             "W":'[ATat]', "w":'[ATat]',
             "K":'[GTgt]', "k":'[GTgt]',
             "M":'[ACac]', "m":'[ACac]',
             "B":'[CGTcgt]', "b":'[CGTcgt]',
             "D":'[AGTagt]', "d":'[AGTagt]',
             "H":'[ACTact]', "h":'[ACTact]',
             "V":'[ACGacg]', "v":'[ACGacg]',
             "N":'[A-Za-z]', "n":'[A-Za-z]'}
                      
    
    with open(motifs_file, 'r') as motifs:
        search_terms = [] 
        un_terms = []
        line = motifs.readline()  
        
        while line:  #each motif
            st = ''
            mt = str(line).strip('\n')
            un_terms.append(mt)
            
            for char in mt:  #each character
                st = st+iupac[char]  #find in iupac, add it to the current search term
        
            search_terms.append(st)
            line = motifs.readline()
    
    return search_terms, un_terms
    
def Trimmer(fasta):
    '''Takes in gene sequence and trims the intronic sequences surrounding the exon sequence to the correct window size.'''
    
    window_size = 300
    exon_bounds = fasta.exon_bounds()
    
    if len(fasta.seq[:exon_bounds[0]]) > window_size: 
        search_start = len(fasta.seq[:exon_bounds[0]])-window_size
    else:
        search_start = 0
    
    if len(fasta.seq[exon_bounds[1]:]) > window_size:
        search_end = exon_bounds[1]+window_size
    else:
        search_end = fasta.tot_seq_len()
    
    trimmed_seq = fasta.seq[search_start:search_end]
    
    return trimmed_seq
        
def Searcher(fasta, motifs):
    '''Takes in gene sequence and motifs, then searches the output of Trimmer and returns the motif positions stored in a dictionary.'''
    
    w_seq = Trimmer(fasta) #trim sequence
    
    motif_positions = {} #find the motifs in the trimmed seq
    for i in range(0,len(motifs[0])):
        motif_positions[motifs[0][i]]=[x.start(0) for x in re.finditer(motifs[0][i], w_seq)]
        motif_positions[motifs[0][i]].append(motifs[1][i])

    return motif_positions
    
def draw(search_results, Motifs):
    '''Draw an image with many genes and many motifs per gene. Output is an .svg'''
    
    ngraphs = len(search_results)
    max_size = 0
    for res in search_results:
        if len(Trimmer(res[0])) > max_size:
            max_size = len(Trimmer(res[0]))  #surface width=max_size+margin
    

    rand_color = randomcolor.RandomColor()
    colors = rand_color.generate(hue="random",luminosity="bright", count=len(search_results)*5)
    motif_color = {}
    i=0
    for motif in Motifs[0]:
        motif_color[motif] = webcolors.hex_to_rgb(colors[i])
        i+=1
    
    surface = cairo.SVGSurface("./motifs_marked.svg", max_size+45, (ngraphs*60)+100)
    graph = cairo.Context(surface)
        
    
    for i in range(0,len(search_results)): #for each gene

        graph.set_source_rgb(0, 0, 0) #draw gene line
        graph.set_line_width(1)
        graph.move_to(45,30*(i+1)+35)
        graph.line_to(len(Trimmer(search_results[i][0]))+45,30*(i+1)+35)
        graph.stroke()
        
        graph.set_font_size(12) #draw gene name
        graph.select_font_face("Optima", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        graph.move_to(5,30*(i+1)+4+35)
        graph.show_text(search_results[i][0].gene)
        
        graph.set_line_width(1.5) #draw exon
        exon_coords = search_results[i][0].exon_bounds_rel(Trimmer(search_results[i][0])) 
        graph.rectangle(exon_coords[0]+45,(30*(i+1)+35)-10,exon_coords[1]-exon_coords[0],18)     
        graph.stroke()
        
        graph.set_source_rgb(255, 255, 255) #fill exon
        exon_coords = search_results[i][0].exon_bounds_rel(Trimmer(search_results[i][0])) 
        graph.rectangle(exon_coords[0]+45+.5,(30*(i+1)+35)-10+.5,exon_coords[1]-exon_coords[0]-1,18-1)     
        graph.fill()
        
        col = 1 
        for key in search_results[i][1].keys(): #draw motifs
            graph.set_source_rgb(motif_color[key][0]/255,motif_color[key][1]/255,motif_color[key][2]/255)
            mot_width = key.count('[')
            for loc in search_results[i][1][key]:
                if type(loc) == int:
                    graph.rectangle(loc+45,(30*(i+1)+35)-14,mot_width, 25)
                    graph.fill()
            col = col+1
    
    graph.move_to(15,(ngraphs*50)+20) #draw legend
    graph.set_source_rgb(0, 0, 0)
    graph.set_font_size(10)
    graph.select_font_face("Optima", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    graph.set_font_size(10)
        
    j = 0
    for i in range(0,len(Motifs[0])):
        graph.set_source_rgb(motif_color[Motifs[0][i]][0]/255,motif_color[Motifs[0][i]][1]/255,motif_color[Motifs[0][i]][2]/255)
        
        graph.rectangle(10+((i+j)*15),(ngraphs*50)+15,5,10)
        graph.fill()
        
        graph.move_to(20+((i+j)*15),(ngraphs*50)+23)
        graph.set_source_rgb(0, 0, 0)
        graph.show_text(search_results[0][1][Motifs[0][i]][-1])
        graph.move_to(20+((i+j)*15),(ngraphs*50))
        
        j+=5
        
        
        
    surface.finish()

#######################
#       MAIN          #
#######################

args = get_arguments()

File = open(args.f,'r')
Motifs = iupac_trans(args.m)

seq_results = []
line = File.readline()

while line:
    to_process = []
    to_process.append(line)
    line = File.readline()
    while (line.startswith('>') != True) and line != '':
        to_process.append(line)
        line = File.readline()
    
    results = Searcher(Sequence(to_process), Motifs)
    
    seq_results.append([Sequence(to_process),results])

File.close()
    
# Draw motif graphs:
draw(seq_results, Motifs)
