#!/usr/bin/python3
import pandas as pd
import re

def fastq_paths(file):
    fastq_df = pd.read_csv(file)
    fq1 = fastq_df["fastq"].str.contains(pat="_R1_", regex=True)
