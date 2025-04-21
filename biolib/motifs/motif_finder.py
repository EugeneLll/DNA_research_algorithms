import neigbors as neighbors_lib


def motif_enumeration(dna: list[str], k: int, d: int) -> set[str]:
    patterns = set()
    for i in range(len(dna[0]) - k + 1):
        pattern = dna[0][i : i + k]
        neighbors = neighbors_lib.neighbors(pattern, d)
        for neighbor in neighbors:
            flag = all(map(lambda x: neighbor in x, dna[1:]))
            if flag:
                patterns.add(neighbor)
