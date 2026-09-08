'''
Functions to plot the root length densities for each
days


@author: Gunther Krauss <guntherkrauss@uni-bonn.de>
'''
from matplotlib import pyplot as plt


'''
Initializes the plot

'''
def init_plot(layers=40):
    ylim = 1
    p,a = plt.subplots(1,1)
    a.set_xlim(0,layers)
    a.set_ylim(0,ylim)
    plt.show(False)
    plt.draw()
    #bckg = p.canvas.copy_from_bbox(a.bbox)
    h=a.plot([],[])[0]
    g=a.plot([],[])[0] 
    return (ylim, p, a, g, h)

'''
Plots the root length density for the actual step
'''    
def plot_step(vRLD, rld_s, ylim, p, a, g, h,layers=40):    
    h.set_data(range(0,layers),vRLD.tolist())
    g.set_data(range(0,layers),rld_s)
    ymax = max(vRLD)
         
        #p.canvas.restore_region(bckg)
    a.draw_artist(h)
    a.draw_artist(g)
    #p.canvas.blit(a.bbox)
    if(ymax > ylim or ymax < 0.1*ylim):
        ylim = 1.3*ymax 
        a.set_ylim(0,ylim)   
        a.relim()
        a.autoscale_view()
    p.canvas.draw()
    return (ylim, p, a, g, h)