import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime as dt

df = pd.DataFrame({'task': ['Ruthenate', 'Nickelate', 'Palladate','Thesis'],
                  'start': pd.to_datetime(['01 Jul 2022', '7 Jan 2022', '01 Sep 2024','01 June 2026']),
                  'end': pd.to_datetime(['31 Aug 2023', '31 Aug 2024', '01 June 2026','12 Jul 2027']),
                  'latest': pd.to_datetime(['31 Dec 2023','31 Dec 2025','31 Dec 2026','01 Jan 2028']) ,
                  'completion_frac': [1, 1, 1, 1]})
print(df)

df['days_to_start'] = (df['start'] - df['start'].min()).dt.days
df['days_to_end'] = (df['end'] - df['start'].min()).dt.days
df['days_to_latest'] = (df['latest'] - df['start'].min()).dt.days
df['task_duration'] = df['days_to_end'] - df['days_to_start'] + 1  # to include also the end date
df['task_maximum'] = df['days_to_latest'] - df['days_to_end'] + 1  # to include also the end date
df['completion_days'] = df['completion_frac'] * df['task_duration']
# plt.barh(y=df['task'], width=df['task_duration'], left=df['days_to_start'])
# plt.show()

# 1
fig = plt.figure(figsize=(24,8))
ax = fig.add_axes([0.2,0.2,0.7,0.6])

plt.barh(y=df['task'], width=df['task_duration'], left=df['days_to_start'] + 1, color = 'orange')
plt.barh(y=df['task'], width=df['task_maximum'], left=df['days_to_end'] + 1,hatch="///",  color='none', edgecolor = 'orange')
plt.barh(y=df['task'], width=df['task_duration'], left=df['days_to_start'] + 1, color = 'none', edgecolor = 'black')
plt.barh(y=df['task'], width=df['task_maximum'], left=df['days_to_end'] + 1,  color='none', edgecolor = 'black')
# plt.title('Proposed PhD Timeline - Julian Nickel', fontsize=15)

# 2
plt.gca().invert_yaxis()

# 3
xticks = np.arange(5, df['days_to_end'].max() + 2, 365)

# 4
xticklabels = pd.date_range(start=df['start'].min() + dt.timedelta(days=4), end=df['end'].max()).strftime("%m/%y")
# 5
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels[::365])

# ax.set_xticks(pd.to_datetime(['01 Jan 2022', '01 Jan 2023', '01 Jan 2024','01 Jan 2025','01 Jan 2026','01 Jan 2027']).dt.days)
# ax.set_xticklabels(['01 Jan 2022', '01 Jan 2023', '01 Jan 2024','01 Jan 2025','01 Jan 2026','01 Jan 2027'])
# 6
ax.xaxis.grid(True, alpha=0.5)

# Marking the current date on the chart
ax.axvline(x=500, color='r', linestyle='dashed')
# ax.text(x=500.5, y=11.5, s='17/11', color='r')

ax.spines[['right', 'top','bottom','left']].set_visible(False)
ax.yaxis.set_tick_params(labelsize=16)
ax.xaxis.set_tick_params(labelsize=16)

fig.savefig('../../Plots/Timeline.svg')
plt.show()