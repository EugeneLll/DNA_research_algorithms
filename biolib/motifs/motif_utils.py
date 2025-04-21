def hamming_distance(str1: str, str2: str) -> int:
    str1 = list(str1)
    str2 = list(str2)
    if len(str1) != len(str2):
        assert 0, "Unequal length"
    distance = 0
    for i in len(str1):
        if str1[i] != str2[i]:
            distance += 1
    return distance
