import os,sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from period_analysis_configuration import TOPOLOGY as T, TOPOLOGY_NAME_DECODER as TD

df = pd.read_csv('graphic/timeseries/total-results.csv')

def detect_variation(df):
	return 'NNODES'
col = detect_variation(df)

df = df.sort_values([col])

# Filter df file by only keeping UIDs that  have not been deleted
#filenamesUID = [x.split('.')[0] for x in os.listdir('graphic/timeseries/timeseries-sweep') if 'series.png' in x]
#df = df[df['UID'].isin(filenamesUID)]
print(df[col])

y = df['T'].to_numpy() * df['SAMPLING_FREQ'].to_numpy()
x = df[col].to_numpy()

f, ax = plt.subplots(1,1,figsize=(20,10))
ax.plot(x,y, c='g',lw=5, ls=':')
ax.scatter(x, y, c='b', s=96)
ax.set_xlabel(col)
ax.set_ylabel('Time')
ax.grid()
ax.set_title(f'Synchronization Time vs {col} for a {TD[T]} Graph')
plt.savefig('graphic/image/period_analysis_result.png')
plt.show()

