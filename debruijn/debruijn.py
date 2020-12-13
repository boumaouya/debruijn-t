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
import random
import statistics
import itertools
import networkx as nx
random.seed(9001)

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
    """Prend un seul argument correspondant au fichier fastq et retourne un
     générateur​ de séquences"""
    if (isfile(fastq_file)):
        with open(fastq_file, "r") as fpt:
            for apt in fpt:
                yield next(fpt).strip()
                next(fpt)
                next(fpt)

def cut_kmer(read, kmer_size):
    """Prend une séquence, une taille de k-mer et retourne un ​ générateur de
    k-mer
    """
    for i, chaine in enumerate(read):
        if(len(read)-i >= kmer_size):  
            seq = read[i:i+kmer_size]
            yield seq

def build_kmer_dict(fastq_file, kmer_size):
    """Prend un fichier fastq, une taille k- mer et retourne un
    dictionnaire ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer
    """
    # Choix d'un dictionnaire
    dictionnaire = dict()
    list_ =[]
    list_distinct = []
    # Récupération k-mer d'une séquence
    for apt in read_fastq(fastq_file):
        for bpt in cut_kmer(apt, kmer_size):
            list_.append(bpt)
    # Enlever les redondances
    for item in list_: 
        if item not in list_distinct: 
            list_distinct.append(item) 
    # Compter le nombre d'occurence par k-mer et remplir le dictionnaire
    for i in list_distinct:
        count = 0
        for item1 in list_: 
            if item1==i: 
                count = count + 1
        dictionnaire[i] = count 
    return dictionnaire

def build_graph(kmer_dict):
    """Prend en entrée un dictionnaire de k-mer et créée l’arbre de
    k-mers préfixes et suffixes décrit précédemment
    """
    # Création d'un graphe
    graph = nx.DiGraph()
    # Supprimer les lettres
    for key, value in kmer_dict.items():
        sequen = key
        npt = len(sequen)
        # Supprimer la première lettre
        new_sequen1 = sequen[0:npt-1]
        # Supprimer la dernière lettre
        new_sequen2 = sequen[1:npt]
        graph.add_edge(new_sequen1, new_sequen2, weight = value)
    return graph
  
def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Prend un graphe et une liste de chemin, la variable booléenne
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés et la variable
    booléenne ​ delete_sink_node pour indiquer si les noeuds de sortie seront supprimés et
    retourne un graphe nettoyé des chemins indésirables
    """
    for path in path_list:
        k=0
        for apt in path:
            if apt in graph.nodes():
                if k==0 and delete_entry_node: graph.remove_node(apt)
                if k==len(path)-1 and delete_sink_node: graph.remove_node(apt)
                if k != 0 and k!=len(path)-1:
                    graph.remove_node(apt)
            k = k+1
    return graph

def std(data):
    """Qui prend une liste de valeur, qui retourne l’écart type
    """
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=False, delete_sink_node=False):
    """Prend un graphe, une liste de chemin, une liste donnant la
    longueur de chaque chemin, une liste donnant le poids moyen de chaque chemin,
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés et
    delete_sink_node pour indiquer si les noeuds de sortie seront supprimés et retourne un
    graphe nettoyé des chemins indésirables
    """
    new_list_path=[]
    remove_list=[]
    tap=[]
    # Récupérer les chemins ayant le plus grand weight_avg
    for k,path in enumerate(path_list):
        if weight_avg_list[k] == max(weight_avg_list):
            new_list_path.append(path)
            tap.append(path_length[k])
        else:
            remove_list.append(path)
    # Parmi les chemins ayant le plus grand weight_avg récupérer ceux qui ont le grand chemin
    best_path = []
    for i,list_ in enumerate(new_list_path):
        if tap[i]==max(tap):
            best_path.append(list_)
        else:
            remove_list.append(list_)
    # Choisir d'une façon aléatoire un chemin parmi les meilleurs, s'ils y en a plusieurs
    npt = len(best_path)
    best = best_path[random.randrange(npt)]
    for zpt in best_path:
        if zpt != best:
            remove_list.append(zpt)
    graph1 = remove_paths(graph, remove_list, delete_entry_node, delete_sink_node)
    return graph1
    
def path_average_weight(graph, path):
    """Prend un graphe et un chemin et qui retourne un poids moyen
    """
    list_weight = []
    for i,seq in enumerate(path):
        if i+1<len(path):
            list_weight.append(graph[seq][path[i+1]]['weight'])
    average = statistics.mean(list_weight)
    return average

def solve_bubble(graph, ancestor_node, descendant_node):
    """Prend un graphe, un noeud ancêtre, un nœud descendant et
    retourne un graph nettoyé de la bulle se trouvant entre ces deux nœuds en utilisant les
    fonctions précédemment développées.
    """
    path_list = list(nx.all_simple_paths(graph,ancestor_node,descendant_node))
    path_length = []
    weight_avg_list = []
    # Récupérer la liste des weight_avg pour chaque chemin
    for path in path_list:
        weight = path_average_weight(graph, path)
        weight_avg_list.append(weight)
    for l in path_list:
        path_length.append(len(l))
    
    graph2 = select_best_path(graph, path_list, path_length, weight_avg_list, delete_entry_node=False, delete_sink_node=False)
    return graph2

def simplify_bubbles(graph):
    """Prend un graphe et retourne un graphe sans bulle
    """
    liste= []
    predecesseurs = []
    # Pour chaque noeud trouver ses prédécesseurs et verifier s'il n'y a pas un ancêtre en coummun entre ces prédécesseurs
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) > 1:
            predecesseurs = list(graph.predecessors(node))
            couples = itertools.combinations(predecesseurs,2)
            for couple in couples:
                rep = nx.lowest_common_ancestor(graph,couple[0], couple[1],None)
                if rep is None:
                    continue
                else:
                    liste.append((rep,node))     
        predecesseurs = []
    for l in liste:
        if l[0] in graph.nodes and l[1] in graph.nodes:
            graph = solve_bubble(graph, l[0], l[1])
    return graph

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    """Prend en entrée un graphe et retourne une liste de noeuds
    d’entrée
    """
    # Liste vide pour stocker les noeuds d'entrée
    nodes = []
    # Récupération des noeuds d'entrée
    for node in graph.nodes():
        pred = graph.predecessors(node)
        if len(list(pred))==0:
            nodes.append(node)
    return nodes

def get_sink_nodes(graph):
    """Prend en entrée un graphe et retourne une liste de noeuds de sortie"""
    nodes = []
    for node in graph.nodes():
        pred = graph.successors(node)
        if len(list(pred))==0:
            nodes.append(node)
    return nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """ Prend un graphe, une liste de noeuds d’entrée et une liste de sortie et
    retourne une liste de tuple(contig, taille du contig)
    """
    # Liste vide pour stocker la liste de couple (contig, taille du contig)
    contigs = []
    for apt in starting_nodes:
        for bpt in ending_nodes:
            # Liste des chemins entre a et b
            list_paths = list(nx.all_simple_paths(graph,apt,bpt))

            # Récupérer les contigs et leurs tailles
            for list_ in list_paths:
                for k, lpt in enumerate(list_):
                    if k==0: cont = list_[k]
                    else:
                        cont = cont + lpt[len(lpt)-1]
                contigs.append((cont,len(cont)))
    return contigs

def save_contigs(contigs_list, output_file):
    """Prend une liste de tuple (contig, taille du contig) et un nom de
    fichier de sortie et écrit un fichier de sortie contenant les contigs selon le format fasta
    """
    def fill(text, width=80):
        """Split text with a line return to respect fasta format"""
        return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))
    
    # Ecrire dans le fichier
    with open(output_file, "wt") as f:
        for k,contigs in enumerate(contigs_list):
            f.write(">contig_" + str(k) + " len=" + str(contigs[1]) + "\n" + fill(contigs[0], width=80) + "\n")
            k = k + 1

    
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
    fastq_file = '../data/eva71_two_reads.fq'
    read_fastq(fastq_file)
    #for a in read_fastq(fastq_file): print(a)

    ########## function cut_kmer ###############
    #for a in cut_kmer('AGGYUUKL', 3): print(a)

    ########## function build_kmer_dict ###############
    fastq_file = '../data/eva71_two_reads.fq'
    dictionnaire = build_kmer_dict(fastq_file, 3)

    ########## function build_graph ###############
    graph = build_graph(dictionnaire)

    ########## function get_starting_nodes ###############
    node_entree = get_starting_nodes(graph)

    ########## function get_sink_nodes #################
    node_sortie = get_sink_nodes(graph)

    ########## function get_contigs ##################
    contigs_list = get_contigs(graph, node_entree, node_sortie)
    
    ########## function save_contigs ##################
    output_file = args.output_file
    save_contigs(contigs_list, output_file)

if __name__ == '__main__':
    main()
