import biolib.neigbors as neighbors_lib
import random
from biolib.motifs import motif_utils
from biolib.motifs.motif import Motif


def motif_enumeration(dna: list[str], k: int, d: int) -> set[str]:
    patterns = set()
    for i in range(len(dna[0]) - k + 1):
        pattern = dna[0][i : i + k]
        neighbors = neighbors_lib.neighbors(pattern, d)
        for neighbor in neighbors:
            flag = all(map(lambda x: neighbor in x, dna[1:]))
            if flag:
                patterns.add(neighbor)


def median_string(dna: list[str], k: int) -> tuple[str, int]:
    dist = float("inf")
    median = ""
    for pattern in motif_utils.k_mer_generator(k):
        d = motif_utils.dna_distance(pattern, dna)
        if dist > d:
            dist = d
            median = pattern
    return median, dist


def greedy_motif_search(dna: list[str], k: int) -> Motif:
    t = len(dna)
    best_motifs = Motif([pattern[:k] for pattern in dna])
    for i in range(len(dna[0]) - k + 1):
        motif = dna[0][i : i + k]
        motifs = Motif([motif])
        for j in range(1, t):
            motifi = max(
                [dna[j][i : i + k] for i in range(len(dna[0]) - k + 1)],
                key=lambda x: motifs.prob(x),
            )
            motifs.add_motif(motifi)
        if motifs.score() < best_motifs.score():
            best_motifs = motifs
    return best_motifs


def greedy_motif_search_laplace(dna: list[str], k: int) -> Motif:
    t = len(dna)
    best_motifs = Motif([pattern[:k] for pattern in dna])
    for i in range(len(dna[0]) - k + 1):
        motif = dna[0][i : i + k]
        motifs = Motif([motif])
        for j in range(1, t):
            motifi = max(
                [dna[j][i : i + k] for i in range(len(dna[0]) - k + 1)],
                key=lambda x: motifs.prob_laplace(x),
            )
            motifs = motifs.add_motif(motifi)
        if motifs.score() < best_motifs.score():
            best_motifs = motifs
    return best_motifs


def randomized_motif_search(dna: list[str], k: int) -> Motif:
    rand = [random.randint(0, len(seq) - k) for seq in dna]
    best = [dna[i][rand[i] : rand[i] + k] for i in range(len(rand))]
    best_motif = Motif(best)
    motif = Motif(best.copy())
    while 1:
        motif = motif.motifs(dna)
        if motif.score() < best_motif.score():
            best_motif = motif.copy()
        else:
            return best_motif


def gibbs_sampler(dna: list[str], k: int, n: int) -> Motif:
    rand = [random.randint(0, len(seq) - k) for seq in dna]
    best = [dna[i][rand[i] : rand[i] + k] for i in range(len(rand))]
    best_motif = Motif(best)
    motif = Motif(best.copy())
    for i in range(n):
        j = random.randint(0, len(dna) - 1)
        destr = motif.profile_destr(dna, j)
        try:
            choise = random.choices(range(len(dna[0]) - k + 1), weights=destr)[0]
        except:
            return best_motif
        new_motif = dna[j][choise : choise + k]
        motifs = motif.as_list()
        motifs[j] = new_motif
        motif = Motif(motifs)
        if motif.score() < best_motif.score():
            best_motif = motif.copy()
    return best_motif
