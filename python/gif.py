
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas as pd
nodes_cor = np.column_stack([nodes_organtype, nodes])

plt.rcParams['figure.figsize'] = [12, 12]
path = 'PiafMunch2_PMA1_output.txt'
#pd.read_table(path,sep='\t')
output = pd.read_table(path,sep='\t',header=1)
JS_ST_begin = (len(node_connection)+1)*31
JS_ST_end = JS_ST_begin +(len(node_connection))
output.iloc[100, JS_ST_begin:JS_ST_end]

#JS_ST_1_out=output.iloc[2, JS_ST_begin:JS_ST_end]
#JS_ST_1_out_array = JS_ST_1_out.values

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x =nodes_c_cor[:,0]
y =nodes_c_cor[:,1]
z =nodes_c_cor[:,2]

# Get rid of colored axes planes
# First remove fill
#ax.set_axis_off()
# Get rid of the panes
#ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
#ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

# Get rid of the spines
ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

# Now set color to white (or whatever is "invisible")
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')


ax.set_xticks([]) 
ax.set_yticks([]) 
ax.set_zticks([])

# Bonus: To get rid of the grid as well:
ax.grid(False)

    


for i in range(0,100):
    connectionflow = ax.scatter(x, y, z, s=100, c=output.iloc[i, JS_ST_begin:JS_ST_end] , cmap=cm.coolwarm, alpha=1, vmin=0, vmax=0.005)
    plt.rcParams.update({'font.size': 20})

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_xlim3d(0, 0.03)
    ax.set_ylim3d(0.03 ,0)
    ax.set_zlim3d(-0.03,0.02)



    #ax.elev = 89.9
    #ax.azim = 270.1
    #ax.dist = 8.0
    #ax.view_init(0, 90)

    cbar = fig.colorbar(connectionflow, shrink=0.5, aspect=10 )
    cbar.ax.get_yaxis().labelpad = 25
    cbar.ax.set_ylabel('Water Flow in Sieve Tube (ml / h)', rotation=270)
    filename='step'+str(i)+'.png'
    plt.savefig(filename, dpi=96)


    plt.show()

#fig.savefig("{}.pdf".format('JW_ST'), bbox_inches='tight')


#def animate(i):
#    connectionflow.set_c(output.iloc[i, JS_ST_begin:JS_ST_end])  # update the data
#    return connectionflow


# Init only required for blitting to give a clean slate.
#def init():
#    connectionflow.set_c(output.iloc[1, JS_ST_begin:JS_ST_end])
#    return connectionflow

#ani = animation.FuncAnimation(fig, animate, np.arange(1, 100), init_func=init,
#                              interval=25, blit=True)
#plt.show()
