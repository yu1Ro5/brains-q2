#%%
import pandas as pd
from src import util
df = pd.read_csv("datasets/dataset.csv")
X,y=util.prepare_data(df)
X.columns
# %%
