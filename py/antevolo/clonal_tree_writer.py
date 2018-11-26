import os
import sys
import shutil
import operator

from Bio import SeqIO
from Bio.Seq import Seq

import igraph
from igraph import *

import utils
import dataset
import clonal_tree_utils

################################ VERTEX WRITERS ##################################
# contract: GetColor(v), GetLabel(v), GetTooltip(v)

class CDR3VertexWriter:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = clonal_tree.FullLengthLineage()
        self.dataset = self.full_length_lineage.Dataset()
        self.cdr3s = list(set([self.dataset.GetCDR3BySeqName(seq.id) for seq in self.clonal_tree.SequenceIter()]))

    def _GetVertexCDR3(self, v):
        return self.dataset.GetCDR3BySeqName(self.clonal_tree.GetSequenceByVertex(v).id)

    def GetColor(self, v):
        vertex_cdr3 = self._GetVertexCDR3(v)
        return utils.GetColorByNormalizedValue('jet', float(self.cdr3s.index(vertex_cdr3)) / len(self.cdr3s))

    def GetLabel(self, v):
        vertex_cdr3 = self._GetVertexCDR3(v)
        return str(self.cdr3s.index(vertex_cdr3))

    def GetTooltip(self, v):
        vertex_cdr3 = self._GetVertexCDR3(v)
        return '\"ID:' + str(v) + " CDR3 ID:" + str(self.cdr3s.index(vertex_cdr3)) + '\"'

class SHMDepthVertexWriter:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = clonal_tree.FullLengthLineage()
        self.dataset = self.full_length_lineage.Dataset()
        self._InitLineageSHMs()

    def _InitLineageSHMs(self):
        self.seq_shm_dict = dict()
        self.min_num_shms = sys.maxint
        self.max_num_shms = 0
        for seq in self.clonal_tree.SequenceIter():
            self.seq_shm_dict[seq.id] = 0
            for gene_type in dataset.AnnotatedGene:
                gene_shms = self.dataset.GetSHMsBySeqName(seq.id, gene_type)
                self.seq_shm_dict[seq.id] += len(gene_shms)
            self.max_num_shms = max(self.max_num_shms, self.seq_shm_dict[seq.id])
            self.min_num_shms = min(self.min_num_shms, self.seq_shm_dict[seq.id])

    def GetColor(self, v):
        seq_id = self.clonal_tree.GetSequenceByVertex(v).id
        rel_value = 0
        if self.max_num_shms - self.min_num_shms != 0:
            rel_value = float(self.seq_shm_dict[seq_id] - self.min_num_shms) / (self.max_num_shms - self.min_num_shms)
        return utils.GetColorByNormalizedValue('Greens', rel_value)

    def GetLabel(self, v):
        return str(v)

    def GetTooltip(self, v):
        seq_id = self.clonal_tree.GetSequenceByVertex(v).id
        return '\"ID:' + str(v) + ' ' + str(self.seq_shm_dict[seq_id]) + ' SHM(s)\"'

class AASeqVertexWriter:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = clonal_tree.FullLengthLineage()
        self.aa_dict = dict()
        self.id_aa_dict = dict()
        for seq in self.clonal_tree.SequenceIter():
            aa_seq = str(Seq(seq.seq).translate())
            self.id_aa_dict[seq.id] = aa_seq
            if aa_seq not in self.aa_dict:
                self.aa_dict[aa_seq] = []
            self.aa_dict[aa_seq].append(seq.id)
        self.max_aa_freq = 0
        for aa in self.aa_dict:
            self.max_aa_freq = max(self.max_aa_freq, len(self.aa_dict[aa]))
        self.max_aa_freq -= 1
        self.sorted_aa_list = sorted(self.aa_dict.keys())

    def GetColor(self, v):
        seq_id = self.clonal_tree.GetSequenceByVertex(v).id
        aa_seq = self.id_aa_dict[seq_id]
        norm_value = 0
        if self.max_aa_freq > 0:
            norm_value = float(len(self.aa_dict[aa_seq]) - 1) / self.max_aa_freq
        return utils.GetColorByNormalizedValue('PuRd', norm_value)

    def GetLabel(self, v):
        return str(v)

    def GetTooltip(self, v):
        seq_id = self.clonal_tree.GetSequenceByVertex(v).id
        aa_seq = self.id_aa_dict[seq_id]
        return '\"AA ID: ' +  str(self.sorted_aa_list.index(aa_seq)) + ' MULT: ' + str(len(self.aa_dict[aa_seq])) + '\"'

class UniqueAAColorWriter:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.aa_seqs = set()
        self.id_aa_dict = dict()
        for seq in self.clonal_tree.SequenceIter():
            aa_seq = str(Seq(seq.seq).translate())
            self.aa_seqs.add(aa_seq)
            self.id_aa_dict[seq.id] = aa_seq
        self.aa_seqs = list(self.aa_seqs)      

    def GetColor(self, v):
        vertex_aa = self.id_aa_dict[self.clonal_tree.GetSequenceByVertex(v).id]
        return utils.GetColorByNormalizedValue('jet', float(self.aa_seqs.index(vertex_aa)) / len(self.aa_seqs))

    def GetLabel(self, v):
        return str(v)

    def GetToolTip(self, v):
        vertex_aa = self.id_aa_dict[self.clonal_tree.GetSequenceByVertex(v).id]
        return '\"AA ID: ' + str(self.aa_seqs.index(vertex_aa))

class RepertoireVertexWriter:
    def __init__(self, clonal_tree, repertoire_fasta):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = clonal_tree.FullLengthLineage() 
        self.repertoire_sequences = set()
        for r in SeqIO.parse(repertoire_fasta, 'fasta'):
            self.repertoire_sequences.add(str(r.seq))

    def GetColor(self, v):
        seq = self.clonal_tree.GetSequenceByVertex(v).seq
        if seq in self.repertoire_sequences:
            return 'green'
        return 'white'

    def GetLabel(self, v):
        return str(v)

    def GetTooltip(self, v):
        return str(v)

class MultiplicityVertexWriter:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.dataset = clonal_tree.FullLengthLineage().Dataset()
        self.tree_mults = [self._GetVertexMultiplicity(v) for v in self.clonal_tree.VertexIter()]
        self.min_mult = min(self.tree_mults)
        self.max_mult = max(10, sorted(self.tree_mults)[len(self.tree_mults) - 2])
        if self.max_mult == self.min_mult:
            self.max_mult = max(self.tree_mults)

    def _GetVertexMultiplicity(self, v):
        return self.dataset.GetSeqMultiplicity(self.clonal_tree.GetSequenceByVertex(v).id)

    def GetColor(self, v):
        vertex_mult = min(self.max_mult, self._GetVertexMultiplicity(v))
        return utils.GetColorByNormalizedValue('Greens', float(vertex_mult - self.min_mult) / (self.max_mult - self.min_mult))

    def GetLabel(self, v):
        return '\"' + str(v) + ' mult:' + str(self._GetVertexMultiplicity(v)) + '\"'

    def GetTooltip(self, v):
        return str(v)

################################ EDGE WRITERS ####################################
# contract: GetTooltip(edge), GetWidth(edge)
class SimpleEdgeWriter:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = clonal_tree.FullLengthLineage()

    def GetTooltip(self, e):
        return str(self.clonal_tree.GetWeightByEdge(e))

    def GetWidth(self, e):
        return str(self.clonal_tree.GetWeightByEdge(e))

    def GetColor(self, e):
        return 'black'

class TypeEdgeWriter:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = clonal_tree.FullLengthLineage()
        self.simple_edge_writer = SimpleEdgeWriter(self.clonal_tree)
        self.edge_type_color = {clonal_tree_utils.DirectedEdgeType.DIRECTED : 'black', clonal_tree_utils.DirectedEdgeType.CDR3 : 'blue', clonal_tree_utils.DirectedEdgeType.REVERSE : 'red', clonal_tree_utils.DirectedEdgeType.SIBLING : 'green'}

    def GetTooltip(self, e):
        return self.simple_edge_writer.GetTooltip(e)

    def GetWidth(self, e):
        return self.simple_edge_writer.GetWidth(e)

    def GetColor(self, e):
        stats = self.clonal_tree.GetStatsByEdge(e)
        return self.edge_type_color[stats.edge_type]

############################# VERSATILE CLONAL TREE WRITER #######################
class ClonalTreeWriter:
    def __init__(self, clonal_tree, vertex_writer, edge_writer):
        self.clonal_tree = clonal_tree
        self.vertex_writer = vertex_writer
        self.edge_writer = edge_writer

    def Output(self, output_base):
        output_dot = output_base + '.dot'
        fh = open(output_dot, 'w')
        fh.write('digraph{\n')
        fh.write('graph [overlaps = false]\n')
        # writing vertices
        for v in self.clonal_tree.VertexIter():
            fh.write(str(v) + ' [style = filled, fillcolor = \"' + self.vertex_writer.GetColor(v) + '\", label = ' + self.vertex_writer.GetLabel(v) + ', tooltip = ' + self.vertex_writer.GetTooltip(v) + ']\n')
        # writing edges
        for e in self.clonal_tree.EdgeIter():
            fh.write(str(e[0]) + ' -> ' + str(e[1]) + ' [color = ' + self.edge_writer.GetColor(e) + ', tooltip = ' + self.edge_writer.GetTooltip(e) + ', penwidth = ' + str(self.edge_writer.GetWidth(e)) + ']\n')
        fh.write('}\n')
        fh.close()
        os.system('fdp -Tsvg ' + output_dot + ' -o ' + output_base + '.svg')
        os.system('fdp -Tpdf ' + output_dot + ' -o ' + output_base + '.pdf')
