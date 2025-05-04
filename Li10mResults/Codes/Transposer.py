import pandas as pd


file1 = input()
file2 = input()

with open(file1, "r") as temp_f:
    # get No of columns in each line
    col_count = [len(l.split(",")) for l in temp_f.readlines()]

### Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
column_names = [i for i in range(0, max(col_count))]

### Read csv
df = pd.read_csv(file1, header=None, delimiter=",", names=column_names)

df.T.to_csv(file2, header=None, index=None, sep=",", mode="w")
