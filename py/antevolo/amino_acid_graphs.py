import os
import sys

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import seaborn as sns
sns.set()

import utils
import dataset
import clonal_tree_constructor
import amino_acid_utils
import vj_annotator
import clonal_tree_writer
import mst_algorithms

def CheckAminAcidGraphEdges(clonal_tree, aa_nucl_dist_map, used_aa, aa_dict):
    print "== Checking amino acid graph..."
    for aa in used_aa:
        nucl_ids = aa_dict.GetIdsByAA(aa) 
        print aa + ", mult: " + str(len(nucl_ids))
        nucl_mults = []
        for nid in nucl_ids:
            nucl_mults.append(clonal_tree.Dataset().GetSeqMultiplicity(nid))
        nucl_mults = sorted(nucl_mults, reverse = True)
        print ','.join([str(m) for m in nucl_mults])
    return
    for e in aa_nucl_dist_map:
        reverse_e = (e[1], e[0])
        if reverse_e not in aa_nucl_dist_map:
            continue
        print "== Processing " + str(e) + ' & ' + str(reverse_e)
        print "Amino acid distance: " + str(utils.HammingDistance(used_aa[e[0]], used_aa[e[1]]))
        old_direct_edges = aa_nucl_dist_map[e]
        old_reverse_edges = aa_nucl_dist_map[reverse_e]
        direct_src_nucl_seqs = [clonal_tree.GetSequenceByVertex(e[0]).seq for e in old_direct_edges]
        reverse_dst_nucl_seqs = [clonal_tree.GetSequenceByVertex(e[1]).seq for e in old_reverse_edges]
        print "Direct distances:"
        for e in old_direct_edges:
            print e, utils.HammingDistance(clonal_tree.GetSequenceByVertex(e[0]).seq, clonal_tree.GetSequenceByVertex(e[1]).seq)
        print "Reverse distances:"
        for e in old_reverse_edges:
            print e, utils.HammingDistance(clonal_tree.GetSequenceByVertex(e[0]).seq, clonal_tree.GetSequenceByVertex(e[1]).seq)
        print "Alternative distances:"
        for e1 in old_direct_edges:
            for e2 in old_reverse_edges:
                print e1[0], e2[1], utils.HammingDistance(clonal_tree.GetSequenceByVertex(e1[0]).seq, clonal_tree.GetSequenceByVertex(e2[1]).seq)

def ComputeAminoAcidGraph(clonal_tree, aa_dict):
    used_aa = []
    aa_edges = dict()
    aa_nucl_dist_map = dict()
    for e in clonal_tree.EdgeIter():
        src_id = clonal_tree.GetSequenceByVertex(e[0]).id
        dst_id = clonal_tree.GetSequenceByVertex(e[1]).id
        src_aa = aa_dict.GetAAById(src_id)
        dst_aa = aa_dict.GetAAById(dst_id)
        if src_aa == dst_aa:
            continue
        if src_aa not in used_aa:
            used_aa.append(src_aa)
        src_index = used_aa.index(src_aa)
        if dst_aa not in used_aa:
            used_aa.append(dst_aa)
        dst_index = used_aa.index(dst_aa)
        aa_edge = (src_index, dst_index)
        if aa_edge not in aa_edges:
            aa_edges[aa_edge] = 0
        aa_edges[aa_edge] += 1
        if aa_edge not in aa_nucl_dist_map:
            aa_nucl_dist_map[aa_edge] = []
        aa_nucl_dist_map[aa_edge].append(e)
    CheckAminAcidGraphEdges(clonal_tree, aa_nucl_dist_map, used_aa, aa_dict)
    return used_aa, aa_edges

class TreeSHMs:
    def __init__(self, full_length_lineage, aa_seqs, aa_dict):
        self.full_length_lineage = full_length_lineage
        self.aa_seqs = aa_seqs
        self.aa_dict = aa_dict
        self.shm_edge_map = dict()
        self.edge_shm_map = dict()
        self._InitializeCDRBounds()

    def _InitializeCDRBounds(self):
        first_seq_id = self.aa_dict.GetIdsByAA(self.aa_seqs[0])[0]
        self.cdr1_bounds = self.full_length_lineage.Dataset().GetCDR1BoundsBySeqName(first_seq_id)
        self.cdr1_bounds = (self.cdr1_bounds[0] / 3, self.cdr1_bounds[1] / 3)
        self.cdr2_bounds = self.full_length_lineage.Dataset().GetCDR2BoundsBySeqName(first_seq_id)
        self.cdr2_bounds = (self.cdr2_bounds[0] / 3, self.cdr2_bounds[1] / 3)
        self.cdr3_bounds = self.full_length_lineage.Dataset().GetCDR3BoundsBySeqName(first_seq_id)
        self.cdr3_bounds = (self.cdr3_bounds[0] / 3, self.cdr3_bounds[1] / 3)

    def AddSHM(self, edge, shm):
        if edge not in self.edge_shm_map:
            self.edge_shm_map[edge] = []
        self.edge_shm_map[edge].append(shm)
        if shm not in self.shm_edge_map:
            self.shm_edge_map[shm] = []
        self.shm_edge_map[shm].append(edge)

    def GetRegionForSHM(self, shm):
        if shm.pos < self.cdr1_bounds[0]:
            return utils.StructuralRegion.FR1
        if shm.pos < self.cdr1_bounds[1]:
            return utils.StructuralRegion.CDR1
        if shm.pos < self.cdr2_bounds[0]:
            return utils.StructuralRegion.FR2
        if shm.pos < self.cdr2_bounds[1]:
            return utils.StructuralRegion.CDR2
        if shm.pos < self.cdr3_bounds[0]:
            return utils.StructuralRegion.FR3
        if shm.pos < self.cdr3_bounds[1]:
            return utils.StructuralRegion.CDR3
        return utils.StructuralRegion.FR4

    def SHMIsRepetitive(self, shm):
        return self.GetSHMMultiplicity(shm) > 1

    def EdgeContainsRepetitiveSHMs(self, edge):
        for shm in self.edge_shm_map[edge]:
            if self.SHMIsRepetitive(shm):
                return True
        return False

    def SHMIter(self):
        for shm in sorted(self.shm_edge_map, key = lambda x : x.pos):
            yield shm

    def GetEdgesBySHM(self, shm):
        return self.shm_edge_map[shm] 

    def GetSHMMultiplicity(self, shm):
        return len(set([e[0] for e in self.shm_edge_map[shm]]))

    def GetNumSHMsOnEdge(self, edge):
        return len(self.edge_shm_map[edge])

    def SHMIsNeutral(self, shm):
        reverse_shm = dataset.SHM(shm.pos, shm.dst_n, shm.src_n)
        return reverse_shm in self.shm_edge_map

    def SHMHasReverse(self, shm):
        reverse_shm = dataset.SHM(shm.pos, shm.dst_n, shm.src_n)
        return self.ContainsSHM(reverse_shm)

    def GetReverseSHM(self, shm):
        return dataset.SHM(shm.pos, shm.dst_n, shm.src_n)

    def ContainsSHM(self, shm):
        return shm in self.shm_edge_map

    def Print(self):
        for e in self.edge_shm_map:
            print e, [(str(shm) + ' : ' + str(len(self.shm_edge_map[shm]))) for shm in self.edge_shm_map[e]]

def ComputeAminoAcidSHMs(full_length_lineage, amino_acids, aa_dict, aa_edges):
    tree_shms = TreeSHMs(full_length_lineage, amino_acids, aa_dict)
    for e in aa_edges:
        for i in range(len(amino_acids[e[0]])):
            if amino_acids[e[0]][i] == amino_acids[e[1]][i]:
                continue
            shm = dataset.SHM(i, amino_acids[e[0]][i], amino_acids[e[1]][i])
            tree_shms.AddSHM(e, shm)
    for shm in tree_shms.SHMIter():
        if tree_shms.SHMIsRepetitive(shm):
            print shm, tree_shms.GetSHMMultiplicity(shm)
    return tree_shms

def GetColorByEdge(edge, tree_shms):
    if tree_shms.EdgeContainsRepetitiveSHMs(edge):
        return 'red'
    return 'black'
    
def OutputAminoAcidGraph(amino_acids, aa_edges, aa_dict, tree_shms, output_base):
    dot_fname = output_base + '.dot'
    fh = open(dot_fname, 'w')
    fh.write('digraph{\n')
    for i in range(len(amino_acids)):
        fh.write(str(i) + ' [label = \"AA ID: ' + str(i) + ' mult: ' + str(aa_dict.GetAAMultiplicity(amino_acids[i])) + '\"]\n')
    for e in aa_edges:
        fh.write(str(e[0]) + ' -> ' + str(e[1]) + ' [label = ' + str(tree_shms.GetNumSHMsOnEdge(e)) + ']\n')
    fh.write('}')
    fh.close()
    os.system('dot -Tpdf ' + dot_fname + ' -o ' + output_base + '.pdf')

def OutputGraphBySHM(amino_acids, aa_edges, aa_dict, tree_shms, shm, output_base):
    dot_fname = output_base + '.dot'
    fh = open(dot_fname, 'w')
    fh.write('digraph{\n')
    for i in range(len(amino_acids)):
        fh.write(str(i) + ' [label = \"AA ID: ' + str(i) + ' mult: ' + str(aa_dict.GetAAMultiplicity(amino_acids[i])) + '\"]\n')
    shm_edges = tree_shms.GetEdgesBySHM(shm)
    reverse_shm = tree_shms.GetReverseSHM(shm)
    reverse_shm_edges = []
    if tree_shms.ContainsSHM(reverse_shm):
        reverse_shm_edges = tree_shms.GetEdgesBySHM(reverse_shm)
    for e in aa_edges:
        color = 'black'
        width = 1
        if e in shm_edges:
            color = 'red'
            width = 5
        if e in reverse_shm_edges:
            color = 'blue'
            width = 5
        fh.write(str(e[0]) + ' -> ' + str(e[1]) + ' [label = ' + str(tree_shms.GetNumSHMsOnEdge(e)) + ', color = \"' + color + '\", penwidth = ' + str(width) + ']\n')
    fh.write('}')
    fh.close()
    os.system('dot -Tpdf ' + dot_fname + ' -o ' + output_base + '.pdf')    

def ColorGraphByAbundantAASHMs(amino_acids, aa_edges, aa_dict, tree_shms, output_base):
    repetitive_shms = [shm for shm in tree_shms.SHMIter() if tree_shms.SHMIsRepetitive(shm)]
    for i in range(len(repetitive_shms)):
        shm_suffix = str(repetitive_shms[i].pos) + '-' + repetitive_shms[i].src_n + '-' + repetitive_shms[i].dst_n
        fname = output_base + '_' + shm_suffix
#        shm_color = utils.GetColorByNormalizedValue('jet', float(i) / len(repetitive_shms))
        OutputGraphBySHM(amino_acids, aa_edges, aa_dict, tree_shms, repetitive_shms[i], fname)

def OutputSHMPlot(full_length_lineage, amino_acids, aa_dict, tree_shms, output_fname):
    aa_len = len(amino_acids[0])
    freq = [0] * aa_len
    for shm in tree_shms.SHMIter():
        freq[shm.pos] += tree_shms.GetSHMMultiplicity(shm)
    seq_id = aa_dict.GetIdsByAA(amino_acids[0])[0]
    fig, ax = plt.subplots(1)
    cdr1_bounds = full_length_lineage.Dataset().GetCDR1BoundsBySeqName(seq_id)
    cdr1_bounds = (cdr1_bounds[0] / 3, cdr1_bounds[1] / 3)
    cdr2_bounds = full_length_lineage.Dataset().GetCDR2BoundsBySeqName(seq_id)
    cdr2_bounds = (cdr2_bounds[0] / 3, cdr2_bounds[1] / 3)
    cdr3_bounds = full_length_lineage.Dataset().GetCDR3BoundsBySeqName(seq_id)
    cdr3_bounds = (cdr3_bounds[0] / 3, cdr3_bounds[1] / 3)
    rect_1 = patches.Rectangle((cdr1_bounds[0], 0), cdr1_bounds[1] - cdr1_bounds[0], max(freq), facecolor = 'red', alpha = 0.25)
    ax.add_patch(rect_1)
    rect_2 = patches.Rectangle((cdr2_bounds[0], 0), cdr2_bounds[1] - cdr2_bounds[0], max(freq), facecolor = 'red', alpha = 0.25)
    ax.add_patch(rect_2)
    rect_3 = patches.Rectangle((cdr3_bounds[0], 0), cdr3_bounds[1] - cdr3_bounds[0], max(freq), facecolor = 'red', alpha = 0.25)
    ax.add_patch(rect_3)
    plt.bar(range(aa_len), freq)
    utils.OutputPlotToPdf(output_fname)

def OutputSHMsToTxt(tree_shms, vj_annotator, output_fname):
    fh = open(output_fname, 'w')
    v_gene = vj_annotator.GetAbundantGene(dataset.AnnotatedGene.V)
    j_gene = vj_annotator.GetAbundantGene(dataset.AnnotatedGene.J)
    fh.write('Position\tSrc_AA\tDst_AA\tEdges\tMultiplicity\tRegion\tHas_reverse\tV_gene\tJ_gene\n')
    for shm in tree_shms.SHMIter():
        edge_str = ','.join([str(e[0]) + '-' + str(e[1]) for e in tree_shms.GetEdgesBySHM(shm)])
        fh.write(str(shm.pos) + '\t' + shm.src_n + '\t' + shm.dst_n + '\t' + edge_str + '\t' + str(tree_shms.GetSHMMultiplicity(shm)) + '\t' + tree_shms.GetRegionForSHM(shm).name + '\t' + str(tree_shms.SHMHasReverse(shm)) + '\t' + v_gene + '\t' + j_gene + '\n')
    fh.close()

############################################################################################
def OutputClonalTree(clonal_tree, full_length_lineage, output_base):
    vertex_writer = clonal_tree_writer.MultiplicityVertexWriter(clonal_tree)
    edge_writer = clonal_tree_writer.TypeEdgeWriter(clonal_tree)
    tree_writer = clonal_tree_writer.ClonalTreeWriter(clonal_tree, vertex_writer, edge_writer)
    tree_writer.Output(output_base)

def OutputAAGraphsForAbundantAAs(full_length_lineage, aa_dict, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    seq_iterator = clonal_tree_constructor.AbundantAASequenceIterator(full_length_lineage, 0.0005, 10)
    empty_filter = clonal_tree_constructor.TrivialFilter()
    edge_computer = clonal_tree_constructor.NaiveEdgeComputer()
    tree_computer = mst_algorithms.VertexMultMSTFinder(full_length_lineage) #IGraphMSTFinder()
    tree_constructor = clonal_tree_constructor.ClonalTreeConstructor(full_length_lineage, seq_iterator, empty_filter, edge_computer, tree_computer, 0)
    clonal_trees = tree_constructor.GetClonalTrees()
    print str(len(clonal_trees)) + ' clonal trees were constructed for abundant AA acids in ' + full_length_lineage.id()
    tree_index = 1
    for tree in clonal_trees:
        root_seq_id = tree.RootSeq().id
        tree_aa_mult = aa_dict.GetAAMultiplicity(aa_dict.GetAAById(root_seq_id))
        OutputClonalTree(tree, full_length_lineage, os.path.join(output_dir, "aa" + str(tree_index) + '_mult' + str(tree_aa_mult)))
        tree_index += 1

############################################################################################
def OutputAbundantAAGraphs(full_length_lineages, output_dir, aa_graph_dir):
    for l in full_length_lineages:
        if len(l) < 100:
            continue
#        custom_filter = clonal_tree_constructor.CustomFilter([clonal_tree_constructor.AbundantVJFilter(l), clonal_tree_constructor.AbundantLengthFilter(l)])
        print "== Processing lineage " + l.id() + '...'
        custom_filter = clonal_tree_constructor.CustomFilter([clonal_tree_constructor.AbundantLengthFilter(l)])
        aa_iterator = clonal_tree_constructor.AllAbundantAAsIterator(l, 0.0005, 10)
        edge_computer = clonal_tree_constructor.HGToolEdgeComputer(os.path.join(output_dir, "full_length_lineages"), 'build/release/bin/./ig_swgraph_construct') # TODO: refactor
        tree_computer = mst_algorithms.VertexMultMSTFinder(l) #IGraphMSTFinder()
        tree_constructor = clonal_tree_constructor.ClonalTreeConstructor(l, aa_iterator, custom_filter, edge_computer, tree_computer, 100)
        clonal_trees = tree_constructor.GetClonalTrees()
        if len(clonal_trees) == 0:
            continue
        clonal_tree = clonal_trees[0]
        aa_dict = amino_acid_utils.AminoAcidDict(l)
        used_aa, aa_edges = ComputeAminoAcidGraph(clonal_tree, aa_dict)
        if len(used_aa) == 0 or len(aa_edges) == 0:
            continue
        OutputAAGraphsForAbundantAAs(l, aa_dict, os.path.join(aa_graph_dir, l.id()))
        annotator = vj_annotator.VJGeneAnnotator(clonal_tree)
        tree_shms = ComputeAminoAcidSHMs(l, used_aa, aa_dict, aa_edges)
        tree_shms.Print()
        OutputAminoAcidGraph(used_aa, aa_edges, aa_dict, tree_shms, os.path.join(aa_graph_dir, l.id()))
#        ColorGraphByAbundantAASHMs(used_aa, aa_edges, aa_dict, tree_shms, os.path.join(output_dir, l.id()))
#        OutputSHMPlot(l, used_aa, aa_dict, tree_shms,  os.path.join(output_dir, l.id() + '_aa.pdf'))
        OutputSHMsToTxt(tree_shms, annotator, os.path.join(aa_graph_dir, l.id() + '_shms.txt'))
