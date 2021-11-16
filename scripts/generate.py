import argparse
import random


def generate(n, density, dist=1.0):
    matrix = dict()
    for _ in range(int(float(n * n) * density * 0.5)):
        i = random.randrange(0, n)
        j = random.randrange(0, n)
        v = random.uniform(-dist, dist)
        matrix[(i, j)] = v
        matrix[(j, i)] = v
    return matrix


def save_mtx(path, n, matrix: dict):
    with open(path, "w") as file:
        file.write(f"{n} {n} {len(matrix)}\n")
        for (i, j), v in matrix.items():
            file.write(f"{i + 1} {j + 1} {v}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--size", default=500)
    parser.add_argument("--density", default=0.2)
    parser.add_argument("--dist", default=50.0)
    parser.add_argument("--path", default="a_500.mtx")
    args = parser.parse_args()

    m = generate(args.size, args.density, args.dist)
    save_mtx(args.path, args.size, m)


if __name__ == '__main__':
    main()
