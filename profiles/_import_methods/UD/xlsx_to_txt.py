import pandas as pd
xlsx = pd.read_excel("complist.xlsx")
# Writing to a text file with space-separated values
xlsx.to_csv("complist.txt", sep=" ", index=False, header=True)
