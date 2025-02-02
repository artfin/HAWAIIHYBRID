from random import randint

def random_with_N_digits(n):
    range_start = 10**(n-1)
    range_end = (10**n)-1
    return randint(range_start, range_end)

if __name__ == "__main__":
    s = "int seeds[100] = {"

    M = 100
    for i in range(M):
        s += str(random_with_N_digits(8))
        if i != M - 1: s += ", "

    s += "};"
    print(s)
