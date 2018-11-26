import os
import sys
import operator

from Bio import SeqIO

import igraph 
from igraph import *

import dataset
import amino_acid_utils
import utils
import clonal_tree_utils

######################## SEQUENCE FILTER #############################
class AbundantVJFilter:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self.most_freq_v = self._FindMostAbundantGene(dataset.AnnotatedGene.V)
        self.most_freq_j = self._FindMostAbundantGene(dataset.AnnotatedGene.J)

    def _FindMostAbundantGene(self, gene_type):
        gene_dict = dict() # gene name -> num sequences
        for seq in self.full_length_lineage.FullLengthSeqIdIter():
            gene_name = self.full_length_lineage.Dataset().GetGeneHitBySeqName(seq.id, gene_type)
            base_name = utils.GetBaseGeneName(gene_name)
            if base_name not in gene_dict:
                gene_dict[base_name] = 0
            gene_dict[base_name] += 1
        most_freq_gene = max(gene_dict.iteritems(), key=operator.itemgetter(1))[0]
#        print "Most frequent " + str(gene_type.name) + ' gene: ' + most_freq_gene + ' (' + str(gene_dict[most_freq_gene]) + ' sequences)'
        return max(gene_dict.iteritems(), key=operator.itemgetter(1))[0]

    def SequenceIsGood(self, seq):
        v_gene = self.full_length_lineage.Dataset().GetGeneHitBySeqName(seq.id, dataset.AnnotatedGene.V)
        j_gene = self.full_length_lineage.Dataset().GetGeneHitBySeqName(seq.id, dataset.AnnotatedGene.J)
        return utils.GetBaseGeneName(v_gene) == self.most_freq_v and utils.GetBaseGeneName(j_gene) == self.most_freq_j

class AbundantLengthFilter:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self._InitLengthDict()

    def _InitLengthDict(self):
        self.len_dict = dict()
        for seq in self.full_length_lineage.FullLengthSeqIdIter():
            if len(seq) not in self.len_dict:
                self.len_dict[len(seq)] = 0
            self.len_dict[len(seq)] += 1
        self.most_freq_length = max(self.len_dict.iteritems(), key=operator.itemgetter(1))[0]
#        print "Most frequent sequence length: " + str(self.most_freq_length) + ' (' + str(self.len_dict[self.most_freq_length]) + ' sequences)'

    def SequenceIsGood(self, seq):
        return len(seq) == self.most_freq_length

class TrivialFilter:
    def SequenceIsGood(self, seq):
        return True

class CustomFilter:
    def __init__(self, filters):
        self.filters = filters

    def SequenceIsGood(self, seq):
        for f in self.filters:
            if not f.SequenceIsGood(seq):
                return False
        return True    

#########################  SEQUENCE ITERATOR #########################
class AllSequenceIterator:
    def __init__(self, full_length_lineage):
        self.full_length_lineage = full_length_lineage
        self.seqs.append([seq for seq in FullLengthSeqIdIter])
        
    def __iter__(self):
        for seqs in self.seqs:
            yield seqs

class AbundantAASequenceIterator:
    def __init__(self, full_length_lineage, min_rel_abundance, min_abs_abundance):
        self.full_length_lineage = full_length_lineage
        self.min_rel_abundance = min_rel_abundance # float number from 0 to 1
        self.min_abs_abundance = min_abs_abundance
        self.aa_dict = amino_acid_utils.AminoAcidDict(full_length_lineage)

    def __iter__(self):
        for aa in self.aa_dict:
            #if float(self.aa_dict.GetAAMultiplicity(aa)) / len(self.full_length_lineage) < self.min_rel_abundance or self.aa_dict.GetAAMultiplicity(aa) < self.min_abs_abundance:
            if self.aa_dict.GetAAMultiplicity(aa) < self.min_abs_abundance:
                continue
            seq_ids = self.aa_dict.GetIdsByAA(aa)
            seqs = [self.full_length_lineage.GetFullLengthSequenceByName(seq_id) for seq_id in seq_ids]
            yield seqs

class AllAbundantAAsIterator:
    def __init__(self, full_length_lineage, min_rel_abundance, min_abs_abundance):
        self.full_length_lineage = full_length_lineage
        self.min_rel_abundance = min_rel_abundance # float number from 0 to 1
        self.min_abs_abundance = min_abs_abundance
        self.aa_dict = amino_acid_utils.AminoAcidDict(full_length_lineage)
        abundant_aa_seqs = []
        for aa in self.aa_dict:
            if float(self.aa_dict.GetAAMultiplicity(aa)) / len(full_length_lineage) < self.min_rel_abundance or self.aa_dict.GetAAMultiplicity(aa) < min_abs_abundance:
                continue
#            print 'adding aa with multiplicity ' + str(self.aa_dict.GetAAMultiplicity(aa))
            seq_ids = self.aa_dict.GetIdsByAA(aa)
            seqs = [self.full_length_lineage.GetFullLengthSequenceByName(seq_id) for seq_id in seq_ids]
            abundant_aa_seqs.extend(seqs)
        self.seqs = []
        self.seqs.append(abundant_aa_seqs)

    def __iter__(self):
        for seqs in self.seqs:
            yield seqs

############################# EDGE COMPUTER #############################
class NaiveEdgeComputer:
    def ComputeEdges(self, seqs):
        edge_dict = dict()
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                edge_dict[(i, j)] = utils.HammingDistance(seqs[i].seq, seqs[j].seq)
        return edge_dict

######################## TREE CONSTRUCTOR #############################
class ClonalTreeConstructor:
    def __init__(self, full_length_lineage, seq_iterator, seq_filter, edge_computer, min_tree_size):
        self.full_length_lineage = full_length_lineage
        self.seq_iterator = seq_iterator
        self.seq_filter = seq_filter
        self.edge_computer = edge_computer
        self.min_tree_size = min_tree_size
        self.clonal_trees = []
        self._ProcessLineage()

    def _ProcessLineage(self):
        for seqs in self.seq_iterator:
            filtered_seqs = [seq for seq in seqs if self.seq_filter.SequenceIsGood(seq)]
#            print '  ' + str(len(seqs)) + ' -> ' + str(len(filtered_seqs))
            if len(filtered_seqs) < self.min_tree_size:
                continue
            graph, edge_weights = self._ComputeHammingGraph(filtered_seqs)
            spanning_tree, tree_weights = self._ComputeSpanningTree(graph, edge_weights)
            undirected_tree = clonal_tree_utils.UndirectedClonalTree(self.full_length_lineage, filtered_seqs)
            for e in spanning_tree.get_edgelist():
                undirected_tree.AddEdge(e, tree_weights[e])
            root_finder = clonal_tree_utils.SimpleRootComputer(undirected_tree)
            directed_tree = clonal_tree_utils.DirectedClonalTree(undirected_tree, root_finder.GetRoot())
            self.clonal_trees.append(directed_tree)

    def _ComputeHammingGraph(self, seqs):
        edge_weights = self.edge_computer.ComputeEdges(seqs)
        graph = Graph()
        graph.add_vertices(len(seqs))
        graph.add_edges(edge_weights.keys())
        return graph, edge_weights

    def _ComputeSpanningTree(self, graph, edge_weights):
        spanning_tree = graph.spanning_tree(weights = edge_weights.values(), return_tree = True)
        tree_weights = dict()
        for e in spanning_tree.get_edgelist():
            tree_weights[e] = edge_weights[e]
        return spanning_tree, tree_weights

    def GetClonalTrees(self):
        return self.clonal_trees

