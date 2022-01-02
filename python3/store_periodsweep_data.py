
import os,sys
import pandas as pd



df = pd.read_csv('graphic/timeseries/total-results.csv')

UID = df['UID'].to_numpy()[-1]

try: os.mkdir(f'graphic/timeseries/timeseries-sweep/{UID}')
except: pass
os.system(f'cp graphic/image/maxdistovertime.png graphic/timeseries/timeseries-sweep/{UID}/fit.png')
os.system(f'cp graphic/image/wholeseries.png graphic/timeseries/timeseries-sweep/{UID}.series.png')
os.system(f'cp graphic/timeseries/result.pkl graphic/timeseries/timeseries-sweep/{UID}.result.pkl')
