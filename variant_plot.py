#%matplotlib inline

import sys

import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

print sys.argv
if len(sys.argv) > 1:
    df = pd.read_csv(sys.argv[1])
else:
    df = pd.read_csv("MAPT_ExACScores.csv")

high_frequency = df[df["ALLELE FREQUENCY"]>0.05]
x = list(high_frequency["AA_POS"])
y = list(high_frequency["ALLELE FREQUENCY"])
mutation = list(high_frequency["MUTATION"])

x_dup = [x[0]]
y_dup = [y[0]]
mutation_dup = [mutation[0]]
for aa in range(1, len(x)):
    if x[aa] == x[aa-1]:
        mutation_dup[-1] = mutation_dup[-1] + ', ' + mutation[aa]
    else:
        x_dup.append(x[aa])
        y_dup.append(y[aa])
        mutation_dup.append(mutation[aa])
        
fig = plt.figure()
ax = fig.add_subplot(111)
x = list(df.AA_POS)
y = list(df["ALLELE FREQUENCY"])

plt.plot(x, y)
plt.axhline(y=0.05, color='r')
for i in range(len(x_dup)):
    ax.annotate(mutation_dup[i], xy=[x_dup[i],y_dup[i]], textcoords='data', rotation=70)
plt.grid()
plt.savefig('variant_plot.tiff')
#plt.show()
