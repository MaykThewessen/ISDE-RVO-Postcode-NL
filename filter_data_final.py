import pandas as pd
import os
clear = os.system('cls' if os.name == 'nt' else 'clear')


isde = pd.read_csv("ISDE_Warmtepomp_all_rows.csv", delimiter=',', encoding='utf-8')

print(isde)

isde = isde.drop(isde.columns[[2, 4, 5, 6, 8, 9, 10, 12]], axis=1)

# Convert all columns to integers (if possible)
for col in isde.columns:
    isde[col] = pd.to_numeric(isde[col], errors='ignore').astype('Int64', errors='ignore')




print(isde)


isde.to_csv("ISDE_Warmtepomp_filtered.csv", index=False)
