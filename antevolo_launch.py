import os
import shutil
import sys

sys.path.append('py/antevolo')
import utils
import dataset
import cdr3_clonal_lineage
import full_length_clonal_lineage
import clonal_tree_writer
import amino_acid_utils
import clonal_tree_utils
import clonal_tree_constructor
import amino_acid_graphs

def FilterClonalLineage(cdr3_lineage):
    return cdr3_lineage.NumFullLengthSequences() < 20 #or cdr3_lineage.NumFullLengthSequences() > 2000

def PrepareOutputDirs(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    output_dirs = dict()
#    output_dirs['trees_CDR3_coloring'] = os.path.join(output_dir, 'trees_cdr3')
#    os.mkdir(output_dirs['trees_CDR3_coloring'])
#    output_dirs['trees_SHM_coloring'] = os.path.join(output_dir, 'trees_shm')
#    os.mkdir(output_dirs['trees_SHM_coloring'])
#    output_dirs['trees_AA_coloring'] = os.path.join(output_dir, 'trees_aa')
#    os.mkdir(output_dirs['trees_AA_coloring'])
#    output_dirs['trees_IgReC_coloring'] = os.path.join(output_dir, 'trees_igrec')
#    os.mkdir(output_dirs['trees_IgReC_coloring'])
    output_dirs['amino_acid_trees'] = os.path.join(output_dir, 'amino_acid_trees')
    if os.path.exists(output_dirs['amino_acid_trees']):
        shutil.rmtree(output_dirs['amino_acid_trees'])
    os.mkdir(output_dirs['amino_acid_trees'])
    output_dirs['amino_acid_graphs'] = os.path.join(output_dir, 'amino_acid_graphs')
    if os.path.exists(output_dirs['amino_acid_graphs']):
        shutil.rmtree(output_dirs['amino_acid_graphs'])
    os.mkdir(output_dirs['amino_acid_graphs'])
    return output_dirs

def OutputClonalTree(fl_lineage, output_dirs):
    cdr3_tree_writer = clonal_tree_writer.ClonalTreeWriter(fl_lineage, clonal_tree_writer.CDR3VertexWriter(fl_lineage), clonal_tree_writer.TypeEdgeWriter(fl_lineage))
    cdr3_tree_writer.Output(os.path.join(output_dirs['trees_CDR3_coloring'], fl_lineage.id()))
    shm_tree_writer = clonal_tree_writer.ClonalTreeWriter(fl_lineage, clonal_tree_writer.SHMDepthVertexWriter(fl_lineage), clonal_tree_writer.TypeEdgeWriter(fl_lineage))
    shm_tree_writer.Output(os.path.join(output_dirs['trees_SHM_coloring'], fl_lineage.id()))
    aa_tree_writer = clonal_tree_writer.ClonalTreeWriter(fl_lineage, clonal_tree_writer.AASeqVertexWriter(fl_lineage), clonal_tree_writer.TypeEdgeWriter(fl_lineage))
    aa_tree_writer.Output(os.path.join(output_dirs['trees_AA_coloring'], fl_lineage.id()))
    igrec_tree_writer = clonal_tree_writer.ClonalTreeWriter(fl_lineage, clonal_tree_writer.RepertoireVertexWriter(fl_lineage, repertoire_fasta), clonal_tree_writer.TypeEdgeWriter(fl_lineage))
    igrec_tree_writer.Output(os.path.join(output_dirs['trees_IgReC_coloring'], fl_lineage.id()))

def OutputAATrees(full_length_lineages, output_dir):
    for l in full_length_lineages:
        print "== Processing lineage " + l.id() + '...'
        custom_filter = clonal_tree_constructor.CustomFilter([clonal_tree_constructor.AbundantVJFilter(l), clonal_tree_constructor.AbundantLengthFilter(l)])
        aa_seq_iterator = clonal_tree_constructor.AbundantAASequenceIterator(l, 0.01, 20)
        edge_computer = clonal_tree_constructor.NaiveEdgeComputer()
        tree_constructor = clonal_tree_constructor.ClonalTreeConstructor(l, aa_seq_iterator, custom_filter, edge_computer, 100)
        clonal_trees = tree_constructor.GetClonalTrees()
        for i in range(len(clonal_trees)):
            tree = clonal_trees[i]
            tree_writer = clonal_tree_writer.ClonalTreeWriter(tree, clonal_tree_writer.SHMDepthVertexWriter(tree), clonal_tree_writer.TypeEdgeWriter(tree))
            tree_writer.Output(os.path.join(output_dir, l.id() + '_tree' + str(i) + '_size' + str(tree.NumVertices())))

def OutputAbundantAAGraphs(full_length_lineages, output_dir):
    for l in full_length_lineages:
        print "== Processing lineage " + l.id() + '...'
        custom_filter = clonal_tree_constructor.CustomFilter([clonal_tree_constructor.AbundantVJFilter(l), clonal_tree_constructor.AbundantLengthFilter(l)])
        aa_iterator = clonal_tree_constructor.AllAbundantAAsIterator(l, 0.0005, 10)
        edge_computer = clonal_tree_constructor.NaiveEdgeComputer()
        tree_constructor = clonal_tree_constructor.ClonalTreeConstructor(l, aa_iterator, custom_filter, edge_computer, 100)
        clonal_trees = tree_constructor.GetClonalTrees()
        if len(clonal_trees) == 0:
            continue
        clonal_tree = clonal_trees[0]
        aa_dict = amino_acid_utils.AminoAcidDict(l)
        used_aa = []
        aa_edges = dict()
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
        if len(used_aa) == 0 or len(aa_edges) == 0:
            continue
        dot_fname = os.path.join(output_dir, l.id() + '.dot')
        fh = open(dot_fname, 'w')
        fh.write('digraph{\n')
        for i in range(len(used_aa)):
            fh.write(str(i) + ' [label = \"AA ID: ' + str(i) + ' mult: ' + str(aa_dict.GetAAMultiplicity(used_aa[i])) + '\"]\n')
        for e in aa_edges:
            fh.write(str(e[0]) + ' -> ' + str(e[1]) + ' [label = ' + str(aa_edges[e]) + ']\n')
            shm_str = []
            for i in range(len(used_aa[e[0]])):
                if used_aa[e[0]][i] != used_aa[e[1]][i]:
                    shm_str.append(str(i) + ':' + used_aa[e[0]][i] + '>' + used_aa[e[1]][i])
            print str(e) + ': ' + ','.join(shm_str)
        fh.write('}')
        fh.close()
        os.system('fdp -Tpdf ' + dot_fname + ' -o ' + os.path.join(output_dir, l.id() + '.pdf'))
                
def main(argv):
    divan_test_dir = argv[1]
    output_dir = argv[2]
#    repertoire_fasta = argv[3]

    output_dirs = PrepareOutputDirs(output_dir)

    dataset_info = dataset.DatasetInfo('test1', 'time_point1', 1, 'cell_type1')
    dataset_obj = dataset.Dataset(dataset_info, divan_test_dir)

    cdr3_lineage_constructor = cdr3_clonal_lineage.CDR3LineageConstructor(dataset_obj, output_dir)
    cdr3_lineages = cdr3_lineage_constructor.Construct()

    full_length_lineages = []
    for l in cdr3_lineages:
        if FilterClonalLineage(l):
            continue
#        print '===='
#        print l.id(), len(l), l.NumFullLengthSequences()
        full_length_lineages.append(full_length_clonal_lineage.FullLengthClonalLineage(l))
    print str(len(full_length_lineages)) + " full-length lineages were constructed"

#    OutputAATrees(full_length_lineages, output_dirs['amino_acid_trees'])
    amino_acid_graphs.OutputAbundantAAGraphs(full_length_lineages, output_dirs['amino_acid_graphs'])

if __name__ == '__main__':
    main(sys.argv)
