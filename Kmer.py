import random
import json
import time

nu_mol_opo = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
nu_mol = ('A', 'T', 'G', 'C')


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def random_seq(p_len):
    return ''.join(random.choices(nu_mol, k=p_len))


class index(object):
    def __init__(self, texta, kmer):
        self.txt = texta
        self.kmer = kmer

    def create_json(self):
        qu = input('Do you really want to update (Y/N):\n')
        print('Go !!!')
        self.start = time.perf_counter()
        if qu == 'Y' or qu == 'y':
            self.p_kmer = {}
            q = [(self.txt[i: i + self.kmer], i) for i in range(len(self.txt) - self.kmer + 1) if
                 self.txt[i: i + self.kmer] != 'N' * self.kmer]
            for w in q:
                self.p_kmer[w[0]] = self.p_kmer.get(w[0], []) + [w[1]]
            data = self.p_kmer
            with open("SEQ.json", "w", encoding="utf8") as write_file:
                json.dump(data, write_file)
                print('JSON DONE !')

    def open_json(self):
        with open("SEQ.json", "r", encoding="utf8") as read_file:
            return json.load(read_file)

    def findPattern(self, pt):
        try:
            return self.json_file[pt]
        except:
            return []

    def findMeSeq(self, nrg):
        self.json_file = self.open_json()
        seq_dict = {}
        counter = 0
        while counter < (4 ** nrg):
            p = random_seq(self.kmer)
            op_p = ''.join(nu_mol_opo[u] for u in p)
            if f'{p}/{op_p}' in seq_dict and counter < 4 ** self.kmer:
                continue
            else:
                seq_dict[f'{p}/{op_p}'] = [self.findPattern(p), self.findPattern(op_p)]
                counter += 1
        return seq_dict

    def printResults(self, rng):
        sequence = self.findMeSeq(rng)
        for key, value in sequence.items():
            if value[0] == [] and value[1] == []:
                continue
            else:
                # print(key, f'{len(value[0]) , len(value[1])}', value)
                if len(value[0]) > 100 or len(value[1]) > 100:
                    print(key, f'{len(value[0]), len(value[1])}', value)

        finish = time.perf_counter()
        print(f'Finished in: {round((finish - self.start), 2)} Seconds')

    def specific(self, ptr):
        start = time.perf_counter()
        op_p = ''.join(nu_mol_opo[u] for u in ptr)
        self.json_file = self.open_json()
        pt_a, pt_b = self.findPattern(ptr), self.findPattern(op_p)
        print(ptr, len(pt_a), pt_a)
        print(op_p, len(pt_b), pt_b)
        finish = time.perf_counter()
        print(f'Finished in: {round((finish - start), 2)} Seconds')


txt = readGenome('seq.txt')
k_mer = 10
#
run = index(txt, k_mer)
run.create_json()
for _ in range(10):
    run.printResults(8)


# txt = readGenome('seq.txt')
# k_mer_seq = 'CTAACCCTAA'
# k_mer = len(k_mer_seq)
# run = index(txt, k_mer)
# run.specific(k_mer_seq)
#
