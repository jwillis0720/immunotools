import os
import sys
import shutil

from Bio import SeqIO
from Bio.Seq import Seq

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt

import dataset
import utils

########################################################################################
class VJGeneAnnotator:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = clonal_tree.FullLengthLineage()
        self.gene_type_mults = dict()
        for gene_type in dataset.AnnotatedGene:
            self.gene_type_mults[gene_type] = dict()
        for gene_type in self.gene_type_mults:
            self._FindMostAbundantGene(gene_type)
            print gene_type, self.gene_type_mults[gene_type]

    def _FindMostAbundantGene(self, gene_type):
        gene_dict = dict() # gene name -> num sequences
        for seq in self.clonal_tree.SequenceIter():
            gene_name = self.full_length_lineage.Dataset().GetGeneHitBySeqName(seq.id, gene_type)
            base_gene = utils.GetBaseGeneName(gene_name)
            if base_gene not in gene_dict:
                gene_dict[base_gene] = 0
            gene_dict[base_gene] += 1
        self.gene_type_mults[gene_type] = gene_dict

    def RealignmentIsNeeded(self):
        for gene_type in self.gene_type_mults:
            if len(self.gene_type_mults[gene_type]) != 1:
                return True
        return False

    def GetAbundantGene(self, gene_type):
        return self.gene_type_mults[gene_type].keys()[0]

class VJUsageAnalyzer:
    def __init__(self):
        self.gene_usage = dict()
        for gene in dataset.AnnotatedGene:
            self.gene_usage[gene] = dict()
        self.lineages = []

    def _UpdateGeneDict(self, gene_type, lineage):
        root_seq_id = lineage.RootSeqId()
        gene_name = utils.GetBaseGeneName(lineage.Dataset().GetGeneHitBySeqName(root_seq_id, gene_type))
        if gene_name not in self.gene_usage[gene_type]:
            self.gene_usage[gene_type][gene_name] = []
        self.gene_usage[gene_type][gene_name].append(lineage.Dataset().GetSHMsBySeqName(root_seq_id, gene_type))

    def AddLineage(self, full_length_lineage):
        self.lineages.append(full_length_lineage)
        for gene in dataset.AnnotatedGene:
            self._UpdateGeneDict(gene, full_length_lineage)

    def Output(self):
        for v in self.gene_usage[dataset.AnnotatedGene.V]:
            print v, [len(shms) for shms in self.gene_usage[dataset.AnnotatedGene.V][v]]

    def _OutputSHMsForV(self, v_gene, output_fname):
        max_length = 300
        num_roots = len(self.gene_usage[dataset.AnnotatedGene.V][v_gene])
        if num_roots < 5:
            return
        for shms in self.gene_usage[dataset.AnnotatedGene.V][v_gene]:
            for shm in shms:
                max_length = max(max_length, shm.pos)
        pos_mult = [0] * max_length
        for shms in self.gene_usage[dataset.AnnotatedGene.V][v_gene]:
            for shm in shms:
                if shm.IsSubstitution():
                    pos_mult[shm.pos] += 1
        plt.bar(range(max_length), pos_mult)
        plt.title(v_gene + ', ' + str(num_roots) + ' roots')
        plt.xlabel('Position (nt)')
        plt.ylabel('# roots')
        utils.OutputPlotToPdf(output_fname)

    def OutputRootSHMs(self, output_dir):
        for v in self.gene_usage[dataset.AnnotatedGene.V]:
            print "Output SHMs for " + v + " gene..."
            self._OutputSHMsForV(v, os.path.join(output_dir, v + '.pdf'))
        
