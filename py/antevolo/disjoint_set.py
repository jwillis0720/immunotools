class DisjointSet:
    def __init__(self, num_elements):
        self._disjoint_set = dict()
        for i in range(num_elements):
            self._disjoint_set[i] = [i]

    def find_index(self, elem):
        for i in self._disjoint_set:
            item = self._disjoint_set[i]
            if elem in item:
                return i
        return None

    def union(self,elem1, elem2):
        index_elem1 = self.find_index(elem1)
        index_elem2 = self.find_index(elem2)
        if index_elem1 != index_elem2 and index_elem1 is not None and index_elem2 is not None:
            self._disjoint_set[index_elem2] = self._disjoint_set[index_elem2] + self._disjoint_set[index_elem1]
            self._disjoint_set.pop(index_elem1)
        return self._disjoint_set
		
    def get(self):
        return self._disjoint_set
