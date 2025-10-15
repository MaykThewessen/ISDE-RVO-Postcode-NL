import pandas as pd


# ---------------------------
# 1. Read ISDE subsidy data
# ---------------------------
# Ensure file path or URL to ISDE Excel file


# Try reading CSV with different delimiters (semicolon is common in European CSVs)
try:
    isde = pd.read_csv("ISDE downloadbestand augustus 2025.csv", delimiter=';', encoding='utf-8')
except:
    # If semicolon doesn't work, try comma or let pandas detect automatically
    try:
        isde = pd.read_csv("ISDE downloadbestand augustus 2025.csv", sep=None, engine='python')
    except:
        # Fall back to reading Excel file directly
        isde = pd.read_excel("ISDE downloadbestand augustus 2025.xlsx")


isde_hp = isde[isde['TECHNIEK'] == 'Warmtepomp']
isde_hp.to_csv("ISDE_Warmtepomp.csv", index=False)

