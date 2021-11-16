import argparse
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).parent.parent
DATA = ROOT / "data"
EXEC = "qr_algo"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exec-folder", default="cmake-build-release")
    args = parser.parse_args()

    threads = [32, 16, 8, 4, 2, 1]
    samples = [20, 5]
    iterations = [10, 15]
    files = ["a_500.mtx", "a_2500.mtx"]

    for f, s, i in zip(files, samples, iterations):
        for t in threads:
            print(f"Exec samples={s} threads={t} iterations={i} file={f}")
            subprocess.check_call([str(ROOT / args.exec_folder / EXEC), str(s), str(t), str(i), str(ROOT / DATA / f)])


if __name__ == '__main__':
    main()
