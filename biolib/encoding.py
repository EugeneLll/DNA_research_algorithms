def pattern_to_number(pattern):
    char_dict = {"a": 0, "c": 1, "g": 2, "t": 3}
    exp = 1
    number = 0
    for i in pattern[::-1]:
        number += exp * char_dict[i]
        exp *= 4
    return number


def number_to_pattern(number, len):
    char_dict = ["a", "c", "g", "t"]
    pattern = ""
    for i in range(len):
        pattern += char_dict[number % 4]
        number = (number - number % 4) // 4
    return pattern[::-1]
