import os
import sys
import Queue

class DFSVertexOrderFinder: 
    def __init__(self, clonal_graph):
        self.clonal_graph = clonal_graph

    def GetOrder(self):
        root = self.clonal_graph.GetRootIndex()
        queue = Queue.PriorityQueue()
        queue.put((0, root))
        order = []
        processed_vertices = set()
        while not queue.empty():
            cur_item = queue.get()
            order.append(cur_item[1])
            processed_vertices.add(cur_item[1])
            child_vertices = self.clonal_graph.GetDescendants(cur_item[1])
            if len(child_vertices) == 1:
                queue.put((cur_item[0], child_vertices[0]))
                continue
            child_prior = cur_item[0]
            for v in child_vertices:
                if v in processed_vertices:
                    continue
                child_prior -= 1
                queue.put((child_prior, v))
        return order

def GetLevelsByVertexOrder(clonal_graph, vertex_order):
    levels = []
    cur_level = 0
    levels.append(cur_level)
    for i in range(0, len(vertex_order) - 1):
        child_vertices = clonal_graph.GetDescendants(vertex_order[i])
        if not(len(child_vertices) == 1 and child_vertices[0] == vertex_order[i + 1]):
            cur_level += 1
        levels.append(cur_level)
    return levels
