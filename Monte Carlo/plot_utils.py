import os
import pickle
from collections import Counter
from datetime import datetime
import plotly.graph_objects as go
import numpy as np


def array_plotter(arr):
    z, x, y = arr.nonzero()
    ele = arr[x, y, z]

    fig = go.Figure(data=[go.Scatter3d(
        x=x,
        y=y,
        z=z,
        opacity=1.0,
        marker=dict(color=ele),
        mode='markers'
    )])
    fig.show()

def retrieve_latest_dump(folder_path="./dump", filter="300"):
    lfolder = os.listdir(folder_path)
    lfolder = [x for x in lfolder if filter in x]
    ans = [datetime.fromisoformat(os.path.splitext(i.split('_')[-1])[0]) for i in lfolder]
    latest = max(ans)
    result = [v for v in lfolder if str(latest) in v]
    print(f"Opening {result[0]}")
    with open(f'./dump/{result[0]}', 'rb') as handle:
        b = pickle.load(handle)

    return b

def sro(mol_frac_dict, key, prob):
    """
                                    alpha_ij = 1 - Pij/cj
    """
    return np.round(1 - prob / mol_frac_dict[key[0]], 3)

def sro_list(arr, neighbour_list, ele_dict):
    sro_counter = {}
    for point in neighbour_list:
        neighbours = [tuple(sorted((arr[point[0]], arr[i]))) for i in point[1]]
        counter = dict(Counter(neighbours))
        total_bonds = sum(list(counter.values()))
        counter = {key: sro(ele_dict, key, value / total_bonds) for key, value in counter.items()}
        for key, value in counter.items():
            if key not in sro_counter:
                sro_counter[key] = [value]
            else:
                sro_counter[key].append(value)

    return sro_counter