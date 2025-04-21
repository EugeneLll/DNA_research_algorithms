import numpy as np


class motif:
    def __init__(pattern_list: list[str]):
        __motifs = np.array(list(map(lambda x: list(x), pattern_list)))
        __dna_parts = ["a", "c", "t", "g"]

    def score(self) -> int:
        pass
