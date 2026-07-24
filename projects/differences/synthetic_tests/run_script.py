#Author: Adam Zelezny
#Date: 2026-07-06

import pyreadr
result = pyreadr.read_r('C:/Users/adamz/Documents/lab/data/GSE231800_S4F_50.rds.gz')
df = result[None]
print(df.info())
print(df.head())
