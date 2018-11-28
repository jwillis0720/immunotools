import sys

class DisjointSet:
    def __init__(self, element_list):
        self._disjoint_set = dict()
        for elem in element_list:
            self._disjoint_set[elem] = set([elem])

    def find_index(self, elem):
        for i in self._disjoint_set:
            item = self._disjoint_set[i]
            if elem in item:
                return i
        print "ERROR: Element " + str(elem) + ' was not found'
        sys.exit(1)
        return None

    def union(self, elem1, elem2):
        index_elem1 = self.find_index(elem1)
        index_elem2 = self.find_index(elem2)
        if index_elem1 != index_elem2:
            for elem in self._disjoint_set[index_elem1]:
                self._disjoint_set[index_elem2].add(elem)
            self._disjoint_set.pop(index_elem1)
		
    def get(self):
        return self._disjoint_set
