class LowFixedAbundanceLeafRemover:
    def __init__(self, clonal_tree, dataset, min_abundance):
        self.clonal_tree = clonal_tree
        self.dataset = clonal_tree.Dataset()
        self.min_abundance = min_abundance

    def VertexToBeRemoved(self, v):
        if not self.clonal_tree.IsRoot(v) and not self.clonal_tree.IsLeaf(v):
            return False
        if self.clonal_tree.IsRoot(v) and len(self.clonal_tree.GetVertexNeighs(v)) > 1:
            return False
        vertex_seq_id = self.clonal_tree.GetSequenceByVertex(v).id
        return self.dataset.GetSeqMultiplicity(vertex_seq_id) < self.min_abundance

class IterativeTipRemover:
    def __init__(self, directed_clonal_tree, vertex_filter):
        self.clonal_tree = directed_clonal_tree
        self.vertex_filter = vertex_filter

    def _CleanTipsInCurrentTree(self):
        vertices_to_remove = []
        for v in self.clonal_tree.VertexIter():
            if self.vertex_filter.VertexToBeRemoved(v):
                vertices_to_remove.append(v)
        for v in vertices_to_remove:
            self.clonal_tree.RemoveVertex(v)
#        print "removed vertices: " + str(vertices_to_remove)
        return len(vertices_to_remove)

    def CleanTips(self):
#        print "Edges: " + str([e for e in self.clonal_tree.EdgeIter()])
        num_removed_vertices = 1
        iteration_step = 1
        while num_removed_vertices != 0: 
            num_removed_vertices = self._CleanTipsInCurrentTree()
            print 'Iteration ' + str(iteration_step) + ': ' + str(num_removed_vertices) + ' vertices were removed'
            iteration_step += 1
        return self.clonal_tree

    
