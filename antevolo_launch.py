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
#    output_dirs['amino_acid_trees'] = os.path.join(output_dir, 'amino_acid_trees')
#    if os.path.exists(output_dirs['amino_acid_trees']):
#        shutil.rmtree(output_dirs['amino_acid_trees'])
#    os.mkdir(output_dirs['amino_acid_trees'])
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

def main(argv):
    divan_test_dir = argv[1]
    output_dir = argv[2]

    output_dirs = PrepareOutputDirs(output_dir)

    dataset_info = dataset.DatasetInfo('test1', 'time_point1', 1, 'cell_type1')
    dataset_obj = dataset.Dataset(dataset_info, divan_test_dir)

    cdr3_lineage_constructor = cdr3_clonal_lineage.CDR3LineageConstructor(dataset_obj, output_dir)
    cdr3_lineages = cdr3_lineage_constructor.Construct()

    full_length_lineages = []
    for l in cdr3_lineages:
        if FilterClonalLineage(l):
            continue
        full_length_lineages.append(full_length_clonal_lineage.FullLengthClonalLineage(l))
    print str(len(full_length_lineages)) + " full-length lineages were constructed"
    amino_acid_graphs.OutputAbundantAAGraphs(full_length_lineages, output_dir, output_dirs['amino_acid_graphs'])

if __name__ == '__main__':
    main(sys.argv)
