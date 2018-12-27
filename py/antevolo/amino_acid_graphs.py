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
import clonal_tree_simplification
import clonal_tree_stats

class ClonalGraph:
    def __init__(self, clonal_tree):
        self.clonal_tree = clonal_tree
        self.full_length_lineage = self.clonal_tree.FullLengthLineage()
        self.vj_annotator = vj_annotator.VJGeneAnnotator(self.clonal_tree)
        self._InitEdges()

    def _InitEdges(self):
        self.aa_dict = amino_acid_utils.AminoAcidDict()
        self.aa_dict.AddClonalTree(clonal_tree)
        self.used_aa = []
        self.aa_edges = dict()
        self.aa_nucl_dist_map = dict()
        for e in self.clonal_tree.EdgeIter():
            src_id = self.clonal_tree.GetSequenceByVertex(e[0]).id
            dst_id = self.clonal_tree.GetSequenceByVertex(e[1]).id
            src_aa = self.aa_dict.GetAAById(src_id)
            dst_aa = self.aa_dict.GetAAById(dst_id)
            if src_aa == dst_aa:
                continue
            src_index = self.aa_dict.GetIndexByAA(src_aa)
            dst_index = self.aa_dict.GetIndexByAA(dst_aa)
            aa_edge = (src_index, dst_index)
            if aa_edge not in self.aa_edges:
                self.aa_edges[aa_edge] = 0
            self.aa_edges[aa_edge] += 1
            if aa_edge not in self.aa_nucl_dist_map:
                self.aa_nucl_dist_map[aa_edge] = []
            self.aa_nucl_dist_map[aa_edge].append(e)

    def ComputeGraphSHMs(self):
        self.shms = TreeSHMs(self.clonal_tree.FullLengthLineage(), self.aa_dict)
        for e in self.aa_edges:
            src_aa = self.aa_dict.GetAAByIndex(e[0])
            dst_aa = self.aa_dict.GetAAByIndex(e[1])
            for i in range(len(src_aa)):
                if src_aa[i] == dst_aa[i]:
                    continue
                shm = dataset.SHM(i, -1, src_aa[i], dst_aa[i])
                self.shms.AddSHM(e, shm)
        for shm in shms.SHMIter():
            if self.shms.SHMIsRepetitive(shm):
                print shm, self.shms.GetSHMMultiplicity(shm)
        return self.shms

    def OutputGraphStructure(self, output_base):
        dot_fname = output_base + '.dot'
        fh = open(dot_fname, 'w')
        fh.write('digraph{\n')
        for aa in self.aa_dict:
            aa_index = self.aa_dict.GetIndexByAA(aa)
            fh.write(str(aa_index) + ' [label = \"AA ID: ' + str(aa_index) + ' mult: ' + str(self.aa_dict.GetAAMultiplicity(aa)) + '\"]\n')
        for e in self.aa_edges:
            fh.write(str(e[0]) + ' -> ' + str(e[1]) + ' [label = ' + str(self.shms.GetNumSHMsOnEdge(e)) + ']\n')
        fh.write('}')
        fh.close()
        os.system('dot -Tpdf ' + dot_fname + ' -o ' + output_base + '.pdf')

    def OutputGraphSHMsToTxt(self, output_fname):
        fh = open(output_fname, 'w')
        v_gene = self.vj_annotator.GetAbundantGene(dataset.AnnotatedGene.V)
        j_gene = self.vj_annotator.GetAbundantGene(dataset.AnnotatedGene.J)
        fh.write('Position\tSrc_AA\tDst_AA\tEdges\tMultiplicity\tRegion\tHas_reverse\tV_gene\tJ_gene\n')
        for shm in self.shms.SHMIter():
            edge_str = ','.join([str(e[0]) + '-' + str(e[1]) for e in self.shms.GetEdgesBySHM(shm)])
            fh.write(str(shm.pos) + '\t' + shm.src_n + '\t' + shm.dst_n + '\t' + edge_str + '\t' + str(self.shms.GetSHMMultiplicity(shm)) + '\t' + self.shms.GetRegionForSHM(shm).name + '\t' + str(self.shms.SHMHasReverse(shm)) + '\t' + v_gene + '\t' + j_gene + '\n')
        fh.close()

    def OutputPutativeAASeqs(self, output_fname):
        fh = open(output_fname, 'w')
        fh.write('Index\tSeq\tDiversity\tNucl_mults\n')
        for aa in self.aa_dict:
            aa_ids = self.aa_dict.GetIdsByAA(aa)
            nucl_mults = sorted([self.full_length_lineage.Dataset().GetSeqMultiplicity(seq_id) for seq_id in aa_ids], reverse = True)
            fh.write(str(self.aa_dict.GetIndexByAA(aa)) + '\t' + aa + '\t' + str(self.aa_dict.GetAAMultiplicity(aa)) + '\t' + ','.join([str(m) for m in nucl_mults]) + '\n')
        fh.close()

class TreeSHMs:
    def __init__(self, full_length_lineage, aa_dict):
        self.full_length_lineage = full_length_lineage
        self.aa_seqs = [aa for aa in aa_dict]
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
        reverse_shm = dataset.SHM(shm.pos, -1, shm.dst_n, shm.src_n)
        return reverse_shm in self.shm_edge_map

    def SHMHasReverse(self, shm):
        reverse_shm = dataset.SHM(shm.pos, -1, shm.dst_n, shm.src_n)
        return self.ContainsSHM(reverse_shm)

    def GetReverseSHM(self, shm):
        return dataset.SHM(shm.pos, -1, shm.dst_n, shm.src_n)

    def ContainsSHM(self, shm):
        return shm in self.shm_edge_map

    def Print(self):
        for e in self.edge_shm_map:
            print e, [(str(shm) + ' : ' + str(len(self.shm_edge_map[shm]))) for shm in self.edge_shm_map[e]]

def OutputClonalTree(clonal_tree, full_length_lineage, output_base, vertex_writer):
    if clonal_tree.NumVertices() > 1000:
        return
    #vertex_writer = clonal_tree_writer.LevelMultiplicityVertexWriter(clonal_tree) #clonal_tree_writer.UniqueAAColorWriter(clonal_tree) #.MultiplicityVertexWriter(clonal_tree)
    edge_writer = clonal_tree_writer.SimpleEdgeWriter(clonal_tree) #TypeEdgeWriter(clonal_tree)
    tree_writer = clonal_tree_writer.ClonalTreeWriter(clonal_tree, vertex_writer, edge_writer)
    tree_writer.Output(output_base)

def OutputLineageAminoAcids(full_length_lineage, aa_dict, output_fname):
    fh = open(output_fname, 'w')
    fh.write('Index\tSeq\tDiversity\tNucl_mults\n')
    for aa in aa_dict:
        aa_ids = aa_dict.GetIdsByAA(aa)
        nucl_mults = sorted([full_length_lineage.Dataset().GetSeqMultiplicity(seq_id) for seq_id in aa_ids], reverse = True)
        fh.write(str(aa_dict.GetIndexByAA(aa)) + '\t' + aa + '\t' + str(aa_dict.GetAAMultiplicity(aa)) + '\t' + ','.join([str(m) for m in nucl_mults]) + '\n')
    fh.close()

############################################################################################
def OutputAbundantAAGraphs(full_length_lineages, output_dir, aa_graph_dir):
    for l in full_length_lineages:
        if len(l) < 100:
            continue
        # clonal tree construction step
        print "== Processing lineage " + l.id() + '...'
        custom_filter = clonal_tree_constructor.CustomFilter([clonal_tree_constructor.AbundantLengthFilter(l)]) #, clonal_tree_constructor.AbundantVJFilter(l)])
        seq_iterator = clonal_tree_constructor.AllSequenceIterator(l) # clonal_tree_constructor.AllAbundantAAsIterator(l, 0.0001, 10)
        edge_computer = clonal_tree_constructor.HGToolEdgeComputer(os.path.join(output_dir, "full_length_lineages"), 'build/release/bin/./ig_swgraph_construct') # TODO: refactor
        tree_computer = mst_algorithms.VertexMultMSTFinder(l) #mst_algorithms.IGraphMSTFinder() # mst_algorithms.VertexMultMSTFinder(l)
        tree_constructor = clonal_tree_constructor.ClonalTreeConstructor(l, seq_iterator, custom_filter, edge_computer, tree_computer, 100)
        clonal_trees = tree_constructor.GetClonalTrees()
        if len(clonal_trees) == 0:
            continue
        clonal_tree = clonal_trees[0]
        # writing clonal tree
        stats_writer = clonal_tree_stats.ClonalTreeStatsWriter(clonal_tree)
        stats_writer.Output(os.path.join(output_dir, l.id() + '_before_simpl.txt'))
        OutputClonalTree(clonal_tree, l, os.path.join(output_dir, l.id() + '_before_simpl'), clonal_tree_writer.LevelMultiplicityVertexWriter(clonal_tree))
        # simplification step
        print "# vertices before simplification: " + str(clonal_tree.NumVertices())
        min_vertex_abundance = 10
        leaf_filter = clonal_tree_simplification.LowFixedAbundanceLeafRemover(clonal_tree, l, min_vertex_abundance)
        leaf_remover = clonal_tree_simplification.IterativeTipRemover(clonal_tree, leaf_filter)
        cleaned_tree = leaf_remover.CleanTips()
        print "# vertices after simplification: " + str(cleaned_tree.NumVertices())
        if cleaned_tree.NumVertices() == 0:
            continue
        # writing clonal tree
        stats_writer = clonal_tree_stats.ClonalTreeStatsWriter(cleaned_tree)
        stats_writer.Output(os.path.join(output_dir, l.id() + '_after_simpl.txt'))
        OutputClonalTree(cleaned_tree, l, os.path.join(output_dir, l.id() + '_after_simpl'), clonal_tree_writer.LevelMultiplicityVertexWriter(cleaned_tree))
        OutputClonalTree(cleaned_tree, l, os.path.join(aa_graph_dir, l.id() + '_nucl_tree'), clonal_tree_writer.UniqueAAColorWriter(cleaned_tree))
        # aa graph
        clonal_graph = ClonalGraph(clonal_tree)
        annotator = vj_annotator.VJGeneAnnotator(clonal_tree)
        graph_shms = clonal_graph.ComputeGraphSHMs()
        graph_shms.Print()
        clonal_graph.OutputGraphStructure(os.path.join(aa_graph_dir, l.id() + '_aa_graph'))
        clonal_graph.OutputGraphSHMsToTxt(os.path.join(aa_graph_dir, l.id() + '_shms.txt'))
        clonal_graph.OutputPutativeAASeqs(os.path.join(aa_graph_dir, l.id() + '_seqs.txt'))
