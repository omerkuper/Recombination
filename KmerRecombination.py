ls = ['A', 'T', 'C', 'G']
matrix_db = [ls[i] + ls[y] for i in range(len(ls)) for y in range(len(ls))]


class recombination:
    def __init__(self, k_mer, inx=None):
        self.k_mer = k_mer // 2
        self.flt_pt = k_mer % 2
        self.data_seq = matrix_db
        self.x = 2 if self.flt_pt == 0 else 1
        self.inx = inx

    def look(self, counter=0):
        if counter <= self.k_mer - self.x:
            if counter < self.k_mer - 1:
                self.data_seq = [outer + inner for outer in self.data_seq for inner in matrix_db][: self.inx]
            else:
                self.data_seq = [outer + inner for outer in self.data_seq for inner in ls][: self.inx]
            counter += 1
            return self.look(counter)
        else:
            return self.data_seq


k_mer = 13
run = recombination(k_mer)
result = run.look()
print(result)
print(len(result))
