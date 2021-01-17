import random as rnd
from pathlib import Path

Path("../sequences/random").mkdir(parents=True, exist_ok=True)

min_len = 1000
max_len = 30001
cantidad = 10
chars = ['A','G','T','C']

for len in range(min_len, max_len, 1000):
    for i in range(1,cantidad+1):
        sequence_list = [chars[rnd.randint(0,3)] for _ in range(len)]
        sequence = ''.join(sequence_list)
        sequence_name = f">seq_{len}_{i}"
        f = open(f"../sequences/random/seq_{len}_{i}.fasta", "w")
        f.write(sequence_name)
        f.write("\n")
        f.write(sequence)
        f.close()
