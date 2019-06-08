import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt('results.dat')

#aprint(np.max(data))
minval=int(np.min(data))
maxval=int(np.max(data))
nbins=maxval-minval+1
print('Number of bins:', nbins)
print()

fig,ax = plt.subplots(1)
plt.hist(data, bins=nbins, range=(minval-0.5,maxval+0.5), rwidth=0.9, align='mid',
         facecolor='g',
         edgecolor='black',
         linewidth=1.2,
         alpha=0.75)
plt.xlabel('Score')
plt.ylabel('Probability')
ax.set_yticklabels([])
plt.yticks([])
plt.draw()
plt.pause(0.01)
#plt.YAxis.remove()
input("<Hit Enter To Close>")
plt.close()
