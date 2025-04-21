def immediate_neighbors(pattern: str) -> set[str]:
    nucl = ["a", "c", "g", "t"]
    neighborhood = set()
    for i in range(len(pattern)):
        symbol = pattern[i]
        for n in nucl:
            if n == symbol:
                continue
            neighbor = list(pattern)
            neighbor[i] = n
            neighborhood.add("".join(neighbor))
    return neighborhood


def neighbors(pattern: str, d: int) -> set[int]:
    neighbourhood = {pattern}
    for i in range(d):
        copy = neighbourhood.copy()
        for pat in copy:
            neighbourhood.update(immediate_neighbors(pat))
    return neighbourhood
