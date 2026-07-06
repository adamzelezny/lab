import pyreadr

# 1. Read the RDS file
result = pyreadr.read_r("C:/Users/adamz/Documents/lab/data/GSE231800_S4F_50.rds.gz")

# 2. Extract the data (read_r returns a dictionary; None is the default key for single objects)
df = result[None]

# 3. Inspect the DataFrame
print(df.info())
print(df.head())
