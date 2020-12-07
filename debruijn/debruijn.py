#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq_file):

    if (isfile(fastq_file)):
        with open(fastq_file, "r") as myf:
            for line in myf:
                yield next(myf).strip()
                next(myf)
                next(myf)

def cut_kmer(read, kmer_size):
    
    for i, chaine in enumerate(read):
        if(len(read)-i >= kmer_size):  
            seq = read[i:i+kmer_size]
            yield seq

def build_kmer_dict(fastq_file, kmer_size):
    dictionnaire = dict()
    list_ =[]
    list_distinct = []

    for a in read_fastq(fastq_file):
        for b in cut_kmer(a, kmer_size):
            list_.append(b)

    for item in list_: 
        if item not in list_distinct: 
            list_distinct.append(item) 

    for i in list_distinct:
        count = 0
        for item1 in list_: 
            if item1==i: 
                count = count + 1
        dictionnaire[i] = count 

    return dictionnaire

def build_graph(kmer_dict):

    graph = nx.DiGraph()
    for key, value in kmer_dict.items():
        sequen = key
        n = len(sequen)
        new_sequen1 = sequen[0:n-1]
        new_sequen2 = sequen[1:n]
        graph.add_edge(new_sequen1, new_sequen2, weight = value)
    return graph
  
def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    pass

def get_sink_nodes(graph):
    pass

def get_contigs(graph, starting_nodes, ending_nodes):
    pass

def save_contigs(contigs_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    ######### fuction read_fastq #############
    #fastq_file = '../data/eva71_two_reads.fq'
    #read_fastq(fastq_file)
    #for a in read_fastq(fastq_file): print(a)

    ########## function cut_kmer ###############
    #for a in cut_kmer('AGGYUUKL', 3): print(a)

    ########## function build_kmer_dict ###############
    #fastq_file = '../data/eva71_two_reads.fq'
    #print(build_kmer_dict(fastq_file, 3))

if __name__ == '__main__':
    main()
