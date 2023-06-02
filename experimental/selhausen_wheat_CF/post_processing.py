import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import pandas as pd

param = 'wheat'

fig = plt.figure()

wheat=pd.read_pickle(param+".pkl")
wheat100=wheat.reshape((100,6,6))
m = np.mean(wheat100, axis =0)
s = np.std(wheat100, axis =0)
z_ = m[:,0]
plt.rcParams["figure.figsize"] = [10, 7]
plt.rcParams["figure.autolayout"] = True

plt.plot(m[:,-3],m[:,0], color = 'r')
plt.fill_betweenx(z_, m[:,-3]+ s[:,-3], m[:,-3] - s[:,-3],
              alpha = 0.3)

for i, j in zip(m[:,-3], m[:,0]):
   plt.text(i, j+0.5, '({}, {})'.format(np.round(i,3), j))

plt.title(param)
plt.xlabel('CF (vRLD/pRLD)')
plt.ylabel('Depth $(cm)$')
plt.rc('font', size=18)
plt.rc('axes', labelsize=18)
plt.xlim(0,4)
fig.savefig(param+'.png', bbox_inches='tight',dpi = 300)
plt.show()

df = pd.DataFrame(m[:,1:], columns = ['pRLD','vRLD','CF(vRLD/pRLD)','RL_image','CF(RL_image/vRLD)'])
df.index = z_
df.to_excel(param+".xlsx")
