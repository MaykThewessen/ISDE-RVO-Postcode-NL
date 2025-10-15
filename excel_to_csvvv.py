import pandas as pd


# ---------------------------
# 1. Read ISDE subsidy data
# ---------------------------
# Ensure file path or URL to ISDE Excel file


isde = pd.read_excel("ISDE downloadbestand augustus 2025.xlsx", engine='openpyxl')

isde_hp = isde[isde['TECHNIEK'] == 'Warmtepomp']
isde_hp.to_csv("ISDE_Warmtepomp_2.csv", index=False)

