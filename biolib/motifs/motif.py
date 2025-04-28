import numpy as np
from collections import Counter


class Motif:
    def __init__(self, pattern_list: list[str]):
        self.__motifs = np.array(
            [np.array(list(x), dtype=np.dtypes.StringDType()) for x in pattern_list]
        )
        self.__dna_parts = ["a", "c", "g", "t"]

    def score(self) -> int:

        t = len(self.__motifs)
        count = sum(
            map(
                lambda x: t - Counter(x).most_common(1)[0][1], self.__motifs.transpose()
            )
        )
        return count

    def count(self) -> np.ndarray[int, int]:

        count = np.array(
            [
                list(map(lambda x: line.get(x, 0), self.__dna_parts))
                for line in map(
                    lambda x: dict(Counter(x).most_common(4)), self.__motifs.transpose()
                )
            ]
        ).transpose()
        return count

    def profile_base(self) -> np.ndarray[int, int]:
        t = len(self.__motifs)
        matr = self.count()
        return matr / t

    def profile_laplace(self) -> np.ndarray[int, int]:
        t = len(self.__motifs)
        matr = self.count() + np.ones_like(self.count())
        return matr / t

    def consensus(self) -> list[str]:
        matr = self.count().argmax(axis=0)
        cons = ""
        for i in matr:
            cons += self.__dna_parts[i]
        return cons

    def prob(self, pattern: str) -> int:
        if len(pattern) != len(self.__motifs[0]):
            raise SyntaxError("Unequal length")
        profile = self.profile_base()
        res = 1
        for i in range(len(pattern)):
            ind = self.__dna_parts.index(pattern[i])
            res *= profile[ind, i]
        return res

    def prob_laplace(self, pattern: str) -> int:
        if len(pattern) != len(self.__motifs[0]):
            raise SyntaxError("Unequal length")
        profile = self.profile_laplace()
        res = 1
        for i in range(len(pattern)):
            ind = self.__dna_parts.index(pattern[i])
            res *= profile[ind, i]
        return res

    def __repr__(self):
        return self.__motifs.__repr__()

    def as_list(self):
        lst = ["".join(i) for i in self.__motifs]
        return lst

    def add_motif(self, motif: str):
        ret = self.as_list()
        ret.append(motif)
        return Motif(ret)

    def motifs(self, dna: list[str]):
        new_motif = []
        k = len(self.__motifs[0])
        for j in range(len(self.__motifs)):
            motifi = max(
                [dna[j][i : i + k] for i in range(len(dna[0]) - k + 1)],
                key=lambda x: self.prob_laplace(x),
            )
            new_motif.append(motifi)
        return Motif(new_motif)

    def profile_destr(self, dna: list[str], line: int) -> list[int]:
        new_motif = self.as_list()
        new_motif.pop(line)
        new_motif = Motif(new_motif)
        k = len(self.count()[0])

        probs = []

        for j in range(len(dna[0]) - k + 1):
            pattern = dna[line][j : j + k]
            res = new_motif.prob_laplace(pattern)
            probs.append(res)

        # print(probs)
        return probs

    def copy(self):
        return Motif(self.as_list())
