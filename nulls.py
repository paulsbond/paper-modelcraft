#!/usr/bin/python3

import pandas as pd


def _main():
    results = pd.read_csv("results.csv")
    todo = results[results.isnull().any(axis=1)]
    print(todo)


if __name__ == "__main__":
    _main()
