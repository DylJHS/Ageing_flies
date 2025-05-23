import pandas as pd
import os

file_directory = os.path.dirname("Results/age_scdgea/ct_specific_MAST/")

full_df = pd.DataFrame()

for file in os.listdir(file_directory):
    if file.endswith(".csv"):
        print(file)
        df = pd.read_csv(os.path.join(file_directory, file))
        if len(df.columns) > 8:
            print("More than 8 columns", '\n',
                  df.columns, '\n')
        full_df = pd.concat([full_df, df], axis=0)


print(full_df.head())
    