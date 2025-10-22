import pandas as pd
import pickle
path = "/Users/albaburgosmondejar/Desktop/shuhui_truth_dm/file_4.pkl"
with open(path, "rb") as f:
    data = pickle.load(f)

print(type(data))
if isinstance(data, dict):
    print(data.keys())
elif isinstance(data, list):
    print(len(data), "items")
else:
    print(data)
df = pd.DataFrame(data)
print(df.head())




