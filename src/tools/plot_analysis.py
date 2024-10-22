import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

import matplotlib as mpl

from collections import defaultdict

# Set font sizes globally
mpl.rcParams['font.size'] = 14  # Default font size for all text
mpl.rcParams['axes.labelsize'] = 16  # Font size for axes labels
mpl.rcParams['axes.titlesize'] = 18  # Font size for axes titles
mpl.rcParams['legend.fontsize'] = 14  # Font size for legends
mpl.rcParams['xtick.labelsize'] = 13  # Font size for x-tick labels
mpl.rcParams['ytick.labelsize'] = 13  # Font size for y-tick labels
mpl.rcParams['figure.titlesize'] = 20  # Font size for figure title

data=pd.read_csv("struct_analysis.csv")

plt.xlabel("supercell size")
plt.ylabel(r"n$th$")



data.sort_values(['struct_vol'] , ascending=[False])


dictionary = defaultdict(list)

for vol, stru in zip(data["struct_vol"], data["first_dep_col_ind"]):
    dictionary[vol].append(stru-1)

size_constant = 5


for xe, ye in dictionary.items():
    xAxis = [xe] * len(ye)

    #square it to amplify the effect, if you do ye.count(num)*size_constant the effect is barely noticeable
    sizes = [ye.count(num) * size_constant for num in ye]
    plt.scatter(xAxis, ye, s=sizes)


plt.savefig("js_vs_vol.pdf")
plt.savefig("js_vs_vol.svg")
plt.show()

