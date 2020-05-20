import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

plt.plot([0,1,2,3,4],[1,2,4,6,8])
plt.savefig('test_plot.pdf',format='pdf')
plt.show()
