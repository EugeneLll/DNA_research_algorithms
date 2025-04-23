import encoding as encode
import neigbors as neighbors_lib


def reverse_compliment(pattern: str) -> str:
    compliment_map = {"t": "a", "g": "c", "a": "t", "c": "g"}
    compliment = "".join(list(map(lambda x: compliment_map[x], list(pattern))))[::-1]
    return compliment


def pattern_position(pattern: str, text: str) -> list[int]:
    pos = 0
    indexes = []
    while text.find(pattern, pos) != -1:
        indexes.append(text.find(pattern, pos))
        pos = indexes[-1]
    return indexes


def compute_frequences(text: str, k: int) -> list[int]:
    freq_array = [0] * pow(4, k)
    for i in range(len(text) - k + 1):
        pattern = text[i : i + k]
        ind = encode.pattern_to_number(pattern)
        freq_array[ind] += 1
    return freq_array


def frequent_words(text: str, k: int) -> tuple[list[str], int]:
    patterns = []
    freqmap = compute_frequences(text, k)
    maxfreq = max(freqmap)
    for i in range(len(freqmap)):
        if freqmap[i] == maxfreq:
            patterns.append(encode.number_to_pattern(i, k))
    return patterns, maxfreq


def find_occurances(text: str, pattern: str) -> int:
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i : i + len(pattern)] == pattern:
            count += 1

    return count


def frequent_words_with_mismatches_compliments(
    text: str, k: int, d: int
) -> tuple[list[str], int]:
    been_to = set()
    patterns = [0] * pow(4, k)
    for i in range(len(text) - k + 1):
        pattern = text[i : i + k]
        compliment = reverse_compliment(pattern)
        ind = encode.pattern_to_number(pattern)
        neighbors_pattern = neighbors_lib.neighbors(pattern=pattern, d=d)
        neighbors_compliment = neighbors_lib.neighbors(compliment, d)
        if been_to.issuperset(neighbors_pattern):
            continue
        neigh = neighbors_pattern.union(neighbors_compliment)
        been_to.update(neigh)
        for neighbor in neigh:
            occurances = find_occurances(text, neighbor)

            patterns[ind] += occurances

    max_patttern = max(patterns)

    final_patterns = []
    for ind in range(len(patterns)):
        if patterns[ind] == max_patttern:
            final_patterns.append(encode.number_to_pattern(ind, k))

    return final_patterns, max_patttern
