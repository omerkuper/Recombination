ls = ['A', 'T', 'C', 'G']
matrix_db = [ls[i] + ls[y] for i in range(len(ls)) for y in range(len(ls))]


def match(data_seq, seq):
    return [outer + inner for outer in data_seq for inner in seq]


class recombination:
    def __init__(self, k_mer):
        self.k_mer = k_mer // 2
        self.flt_pt = k_mer % 2
        self.data_seq = matrix_db
        self.x = 2 if self.flt_pt == 0 else 1

    def look(self, count=0, inx=None):
        inx, counter = inx, count
        if counter <= self.k_mer - self.x:
            if counter < self.k_mer - 1:
                self.data_seq = match(self.data_seq, matrix_db)[: inx]
            else:
                self.data_seq = match(self.data_seq, ls)[: inx]
            counter += 1
            return self.look(counter, inx)
        else:
            return self.data_seq


k_mer = 10
run = recombination(k_mer)
result = run.look()
print(result)
print(len(result))

