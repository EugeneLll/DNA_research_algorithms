import itertools


def hamming_distance(str1: str, str2: str) -> int:
    str1 = list(str1)
    str2 = list(str2)
    if len(str1) != len(str2):
        assert 0, "Unequal length"
    distance = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            distance += 1
    return distance


def dna_distance(pattern: str, dna: list[str]) -> int:
    return sum(
        [
            min(
                [
                    hamming_distance(pattern, line[i : i + len(pattern)])
                    for i in range(len(line) - len(pattern) + 1)
                ]
            )
            for line in dna
        ]
    )


def k_mer_generator(k: int) -> str:
    nucl = ["a", "c", "g", "t"]
    return ["".join(i) for i in itertools.product(nucl, repeat=k)]
