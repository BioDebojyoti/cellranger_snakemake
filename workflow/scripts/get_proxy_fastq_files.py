import pandas as pd
import re
import json
import sys


def get_proxy(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)

    folders_expected = pd.DataFrame.from_dict(data, orient="index")[0].tolist()
    folder_proxies = [re.sub("\/$", ".csv", l) for l in folders_expected]

    for file in folder_proxies:
        with open(file, "w"):
            pass


if __name__ == "__main__":
    json_file = sys.argv[1]
    get_proxy(json_file)
