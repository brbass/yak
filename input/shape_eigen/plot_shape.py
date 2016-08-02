import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt("plot_data.txt")
data = data[np.argsort(data[:,1])]
dist_values = np.unique(data[:,1])

plt.figure()
for i in range(len(dist_values)):
    indices = data[:, 1] == dist_values[i]
    local_data = data[indices]
    local_data = local_data[np.argsort(local_data[:, 0])]
    
    plt.plot(local_data[:, 0], local_data[:, 2], label=str(dist_values[i]))
plt.legend(title="min distance")
plt.xlabel("shape parameter")
plt.ylabel("k eigenvalue")
plt.show()
    

# data_ga = data[indices_ga]
# data_mq = data[indices_mq]
# data_imq = data[indices_imq]

# colors = ['#1b9e77','#d95f02','#7570b3']
# shapes = ['^', 'o', 's']
# # colors = iter(['#66c2a5','#fc8d62','#8da0cb'])

# plt.figure()
# ax = plt.subplot(111)
# plt.plot(data_ga[:,0], data_ga[:,4],
#          linestyle='solid', color=colors[0],
#          marker=shapes[0], markeredgecolor=colors[0], markeredgewidth=1, markerfacecolor='None', 
#          label="GA")
# plt.plot(data_mq[:,0], data_mq[:,4],
#          linestyle='solid', color=colors[1],
#           marker=shapes[1], markeredgecolor=colors[1], markeredgewidth=1, markerfacecolor='None',
#          label="MQ")
# plt.plot(data_imq[:,0], data_imq[:,4],
#          linestyle='solid', color=colors[2],
#          marker=shapes[2], markeredgecolor=colors[2], markeredgewidth=1, markerfacecolor='None',
#          label="IMQ")
# plt.xlabel(r"shape multiplier $k$")
# plt.ylabel(r"$L_2$ relative error in scalar flux $\phi$")
# # ax.set_xticks(major_ticks)
# # ax.set_xticks(minor_ticks, minor=True)
# # ax.set_yticks(y_ticks)
# plt.ylim(0, 0.01)
# plt.xlim(0, 3)
# plt.grid(True)
# plt.legend(fontsize=12)
# plt.tight_layout()
# plt.savefig("figures/shape.pdf")
# plt.show()
# plt.close()
