import sys
import datetime
import matplotlib.pylab as plt
import numpy as np
import timeit
import os
import plotly
import os
import math
import plotly.graph_objects as go
import xml.etree.ElementTree as ET
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from plotly.subplots import make_subplots

sys.path.append("../..")
import plantbox as pb #CPlantBox Python Binding
import matplotlib.pylab as plt
import numpy as np
from scipy import stats
import pandas as pd
import datetime

def visual_plant_sub(plant1,name='plant'):
    nodes_cor = python_nodes(plant1) # use the object name created previously to get its coordinates
    fig= go.Scatter3d(
        x=nodes_cor.T[3],
        y=nodes_cor.T[4],
        z=nodes_cor.T[5],
        mode='markers',
        marker=dict(
            size=3,
            color=nodes_cor.T[1],                # nodes_cor.T[1] is organ type, nodes_cor.T[2] is the connection number of a node 
            colorscale=[[0, "wheat"], #color of the root, change it to "yellow" to see the difference
                    [0.5, "darkgreen"],
                    [1.0, "lightgreen"],],  opacity=0.8
        ), name =name
    )

    return fig
def change_parameter(input_name, output_name, organ_name, subtype , parameter_name, value_type, value):   
    all_parameter = ET.parse("../../modelparameter/{}.xml".format(input_name)) # read the parameter file from xml file
    plant_parameter = all_parameter.getroot() # get the first level of parameters
    original_lmax = plant_parameter.find("./organ[@type='{}'][@subType='{}']/parameter[@name='{}']".format(organ_name, subtype, parameter_name)).get('{}'.format(value_type)) # get function to read the value
    # The '10' in the following line is the value that need to be changed change the value
    plant_parameter.find("./organ[@type='{}'][@subType='{}']/parameter[@name='{}']".format(organ_name, subtype, parameter_name)).set('{}'.format(value_type),'{}'.format(value)) # set function to read the value
    current_lmax = plant_parameter.find("./organ[@type='{}'][@subType='{}']/parameter[@name='{}']".format(organ_name, subtype, parameter_name)).get('{}'.format(value_type)) # get function to read the value
    print('original {} of {} organ with subtype {} is {}, changed to {}'.format(parameter_name, organ_name, subtype ,original_lmax, current_lmax))
    all_parameter.write('../../modelparameter/{}.xml'.format(output_name))   

def visual_plant(plant1):
    subfig = visual_plant_sub(plant1)
    fig = make_subplots(
    rows=1, cols=1,
    specs=[[{'type': 'surface'}]])
    fig.add_trace(subfig) 
    fig.update_layout(scene_aspectmode='data',)
    return fig

def CPlantBox_PiafMunch(name, time, output = "test_output"):
    plant = pb.Plant()
    plant.openXML("/modelparameter/plant/" + name)
    seeds = plant.getOrganRandomParameter(pb.OrganTypes.seed)
    roots = plant.getOrganRandomParameter(pb.OrganTypes.root)
    stems = plant.getOrganRandomParameter(pb.OrganTypes.stem)
    leafs = plant.getOrganRandomParameter(pb.OrganTypes.leaf)
    plant.initialize(True)
    plant.simulate(time)
    plant.write("{}.vtp".format(str(output)))
    dict_all  = convert( plant )
    return dict_all ;

def CPlantBox(name, time, output = "output"): #define a function, in line 20, we can run it in one line of code
    plant = pb.Plant()
    plant.openXML('../../modelparameter/plant/'+name)
    seeds = plant.getOrganRandomParameter(pb.OrganTypes.seed)
    roots = plant.getOrganRandomParameter(pb.OrganTypes.root)
    stems = plant.getOrganRandomParameter(pb.OrganTypes.stem)
    leafs = plant.getOrganRandomParameter(pb.OrganTypes.leaf)
    plant.initialize(True)
    plant.simulate(time)
    plant.write("{}.vtp".format(str(output)))
    return plant;

def CPlantBox_analysis(name, time, output = "output"): #define a function, in line 20, we can run it in one line of code
    plant = pb.Plant()
    plant.openXML("../../modelparameter/plant/" + name)
    plant.initialize(True)
    plant.simulate(time)
    #plant.write("../../results/{}.vtp".format(output),15)
    ana = pb.SegmentAnalyser(plant)
    ana.write("../../results/{}.vtp".format(str(output)))
    ana.write("../../results/{}.txt".format(str(output)))
    return plant;

def python_nodes(plant):
    organtype = np.array([np.array(s) for s in plant.getParameter('organType')])
    lines = list([list(s) for s in plant.getPolylines()])
    segs = np.array([np.array(s) for s in plant.getSegments()])
    nods = np.array([np.array(s) for s in plant.getNodes()])
    node_organtype = np.zeros(len(nods)-1)
    k=0
    for i in range(0,len(organtype)):
        for ii in range(0,len(lines[i])-1):
            node_organtype[k]= organtype[i]
            k=k+1
    nodes = nods/100 # convert from cm to m 
    node_connection_o = segs
    nodes_with_organtype = np.column_stack([node_connection_o, node_organtype]) #make the node has organtype, to know the source sink relation
    node_connection1, node_connection2 = np.split(node_connection_o.T,2) #seperate the 2-number nodes into two list
    node_connection1 = np.row_stack([node_connection1, node_organtype]) #first numbe list
    node_connection2 = np.row_stack([node_connection2, node_organtype]) #second number list
    nodes_organtype = np.column_stack([node_connection1,node_connection2])
    _, indices = np.unique(nodes_organtype.T[:,0], return_index=True) #sort the list to remove duplicates
    nodes_organtype = nodes_organtype.T[indices,:]
    nodes_cor = np.column_stack([nodes_organtype, nodes]) # adding coordinates into the connections
    node_connection = np.copy(node_connection_o)
    unq, unq_idx, unq_cnt = np.unique(node_connection, return_inverse=True, return_counts=True)# check if all the connections are unique
    nodes_organtype = np.column_stack((nodes_organtype,unq_cnt ))
    nodes_organtype.astype(np.int_)
    node_connection.astype(np.int_)
    nodes_cor = np.column_stack([nodes_organtype, nodes])
    return nodes_cor;

# Convert the connected nodes in CPlantBox to PiafMunch arbitary numbers.
def convert( plant ):
    organtype = np.array([np.array(s) for s in plant.getParameter('organType')])
    lines = list([list(s) for s in plant.getPolylines()])
    segs = np.array([np.array(s) for s in plant.getSegments()])
    nods = np.array([np.array(s) for s in plant.getNodes()])
    node_organtype = np.zeros(len(nods)-1)
    k=0
    for i in range(0,len(organtype)):
        for ii in range(0,len(lines[i])-1):
            node_organtype[k]= organtype[i]
            k=k+1
    nodes = nods/100 # convert from cm to m 
    node_connection_o = segs
    nodes_with_organtype = np.column_stack([node_connection_o, node_organtype]) #make the node has organtype, to know the source sink relation
    node_connection1, node_connection2 = np.split(node_connection_o.T,2) #seperate the 2-number nodes into two list
    node_connection1 = np.row_stack([node_connection1, node_organtype]) #first numbe list
    node_connection2 = np.row_stack([node_connection2, node_organtype]) #second number list
    nodes_organtype = np.column_stack([node_connection1,node_connection2])
    _, indices = np.unique(nodes_organtype.T[:,0], return_index=True) #sort the list to remove duplicates
    nodes_organtype = nodes_organtype.T[indices,:]
    nodes_cor = np.column_stack([nodes_organtype, nodes]) # adding coordinates into the connections
    node_connection = np.copy(node_connection_o)
    unq, unq_idx, unq_cnt = np.unique(node_connection, return_inverse=True, return_counts=True)# check if all the connections are unique
    nodes_organtype = np.column_stack((nodes_organtype,unq_cnt ))
    nodes_organtype.astype(np.int_)
    node_connection.astype(np.int_)
    nodes_cor = np.column_stack([nodes_organtype, nodes])
    nodes_cor = np.column_stack([nodes_organtype, nodes]) # adding coordinates into the connections
    node_connection = np.copy(node_connection_o)
    unq, unq_idx, unq_cnt = np.unique(node_connection, return_inverse=True, return_counts=True)# check if all the connections are unique
    nodes_organtype = np.column_stack((nodes_organtype,unq_cnt ))
    nodes_organtype.astype(np.int_)
    node_connection.astype(np.int_)
    nodes_cor = np.column_stack([nodes_organtype, nodes])
    stem_nodes = nodes_organtype[(nodes_organtype[:,1]== 8)|(nodes_organtype[:,1]== 4)][:,0]
    index_stem= list(range(1, len(stem_nodes)+1))
    for i in range(len(stem_nodes)):
        index_stem[i] = np.where( node_connection[:,1] == stem_nodes[i])

    for i in range(len(index_stem)-1):
        node_connection[index_stem[i+1][0][0]]=node_connection[index_stem[i+1][0][0]][::-1]

    #print("node_connection length is{}".format(len(node_connection) ))
    # Write input file for PiafMunch
    # Here we created several functions we used later

    pi= 3.1415926

    k = 1.118e-12 # conductivity (m2)
    n = 1.7e-9 / 3600 # viscosity (MPa.h)
    s = (12.880e-6)**2*pi # surface (m2)
    l = 5e-2 # length (m)

    K1 = (k/n)*(s/l) *10**6

    K_1 = (1 / K1)
    #print(K_1)

    # 1 meter inner phloem resistance

    k = 0.693e-12 # conductivity (m2)
    n = 1.7e-9 / 3600 # viscosity (MPa.h)
    s = (14.6e-6)**2*pi # surface (m2)
    l = 5e-2 # length (m)

    K1 = (k/n)*(s/l) *10**6

    K_4 = (1 / K1)
    #print(K_4)
    # 4 meter inner phloem resistance

    k = 0.879e-12 # conductivity (m2)
    n = 1.7e-9 / 3600 # viscosity (MPa.h)
    s = (9.83e-6)**2*pi # surface (m2)
    l = 5e-2 # length (m)

    K1 = (k/n)*(s/l) *10**6

    K_7 = (1 / K1)
    #print(K_7)

    #7 meter inner phloem resistance

    k = 1.795e-12 # conductivity (m**2)
    n = 1.7e-9 / 3600 # viscosity (MPa.h)
    s = (13.830e-6)**2*pi # surface (m**2)
    l = 5e-2 # length (m)

    K1 = (k/n)*(s/l) *10**6 

    K_1_i = (1 / K1)
    #print(K_1_i)
     
    k = 0.98712e-12 # conductivity (m**2)
    n = 1.7e-9 / 3600 # viscosity (MPa.h)
    s = (19.710e-6)**2*pi # surface (m**2)
    l = 5e-2 # length (m)

    K1 = (k/n)*(s/l) *10**6

    K_4_i = (1 / K1)
    #print(K_4_i)
     
    k = 0.3532e-12 # conductivity (m**2)
    n = 1.7e-9 / 3600 # viscosity (MPa.h)
    s = (10.670e-6)**2*pi # surface (m**2)
    l = 5e-2 # length (m)

    K1 = (k/n)*(s/l) *10**6

    K_7_i = (1 / K1)
    #print(K_7_i)

    R_name = [1,4,7]
    R_ex_number = [K_1,K_4,K_7] #external phloem resistance
    R_in_number = [K_1_i,K_4_i,K_7_i] #internal phloem resistance


    plt.rcParams['figure.figsize'] = [10, 10]

    slope, intercept, r_value, p_value, std_err = stats.linregress(R_name[0:2], R_ex_number[0:2]) #least square regression
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(R_name[1:3], R_ex_number[1:3]) #least square regression
    #plt.plot(R_name[0:2], R_ex_number[0:2], 'o', label='Resistance in external phloem')
    #plt.plot(R_name[1:3], R_ex_number[1:3], 'o', label='Resistance in external phloem', c='black' )
    #plt.plot(R_name[0:2], intercept + slope*R_name[0:2], 'r', label='fitted line1')
    #plt.plot(R_name[1:3], intercept2 + slope2*R_name[1:3], 'r', label='fitted line2')
    #plt.xlabel("Distance to seed (meter)")
    #plt.ylabel("Resistance (MPa h mL$^{-1}$))")
    #plt.legend()
    #plt.show()


    stem_nodes = nodes_organtype[(nodes_organtype[:,1]== 4)][:,0] 
    root_nodes = nodes_organtype[(nodes_organtype[:,1]== 2)][:,0] 
    leaf_nodes = nodes_organtype[(nodes_organtype[:,1]== 8)][:,0] 
    def organ_marker(int):
        if int in leaf_nodes:
        #if int == 21: #int is source 
            return 'o'
        elif int in root_nodes : #int is sink
            return 'x'
        #elif int == 40 :
        elif int in stem_nodes:
            return '*'
        else:
            return 'None'
    
    #node_connection_cor
    nodes_c_cor = np.full((len(node_connection), 3),0.0)
    #charar = np.chararray((3, 3))
    nodes_org = np.full((len(node_connection), 1),0)
    for i in range(0,len(node_connection)):
        nodes_org[i]= nodes_organtype[node_connection[i][0]][1]
    nodes_c_marker = np.full((len(node_connection), 1),'k')
    for i in range(0,len(node_connection)):
        nodes_c_marker[i]= organ_marker(node_connection[i][0])
    #nodes_c_marker[0][0]= 'D'
    nodes_c_marker.astype(str, copy = True)


    for i in range(0,len(node_connection)): #calculate the mid point of two segment nodes
        nodes_c_cor[i] = (nodes[node_connection[i,0]] + nodes[node_connection[i,1]])/2


        
        
    b = np.zeros((nodes_organtype.shape[0],nodes_organtype.shape[1]+1)); b[:,:-1] = nodes_organtype
    node_c_o = np.zeros((node_connection.shape[0],node_connection.shape[1]+5)); node_c_o[:,:-5] = node_connection

    for i in range(0,len(node_connection)):
        node_c_o[i][2:5] = (nodes[node_connection[i,0]] + nodes[node_connection[i,1]])/2
    #    node_c_o[i][6] = nodes_organtype[node_c_o[i][0]][0]
        #node_c_o[i][6] = nodes_organtype[node_c_o[i][0]][1]

        


    r_st_all = np.full((len(node_connection), 1),0.0)
    for i in range(0,len(node_connection)):
        r_st_all[i]= intercept + slope*node_c_o[i,4]


        
    node_c_o = pd.DataFrame({'1st_node':node_connection[:,0],'2nd_node':node_connection[:,1],'x':node_c_o[:,2]
                            ,'y':node_c_o[:,3],'z':node_c_o[:,4], 'organ_type':nodes_org[:,0],'marker':nodes_c_marker[:,0], 'r_st':r_st_all[:,0]})
    #claculate the length from seed of every nodes (root is minus)
    nodes_length = np.zeros(len(nodes_cor))
    for i in range(2,len(nodes_cor)):
        if nodes_cor[i][1]==2: #if the organtype is root
            nodes_length[i] = nodes_length[node_connection_o[(node_connection_o[:,1] == i)][0][0]] - ((nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][0]][3]-nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][1]][3])**2 + (nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][0]][4]-nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][1]][4])**2 + (nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][0]][5]-nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][1]][5])**2)**(0.5) 
        else:
            nodes_length[i] = nodes_length[node_connection_o[(node_connection_o[:,1] == i)][0][0]] + ((nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][0]][3]-nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][1]][3])**2 + (nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][0]][4]-nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][1]][4])**2 + (nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][0]][5]-nodes_cor[node_connection_o[(node_connection_o[:,1] == i)][0][1]][5])**2)**(0.5) 

            
            
    #     print(nodes_length[i])
    #     print(i)
        #print(node_connection[i])
        #the length of the first node of connection i in node_connection
        #node_c_o[i][6] = nodes_organtype[node_c_o[i][0]][0]
        #node_c_o[i][6] = nodes_organtype[node_c_o[i][0]][1]

        #manual check the length of all nodes.
    nodes_cor[node_connection[(node_connection[:,1] == 3)][0][0]][2]
    #len(nodes_length)
    #print(len(nodes_length))
    nodes_length
    #node_connection_o[0,0]=1
    np.c_[ nodes_organtype, nodes_length]


    len(nodes_organtype)
    #node_connection
    len(nodes_cor)
    nodes_cor[:,4]
    nodes_r_st = np.zeros(len(nodes_cor))
    for i in range(0,len(nodes_cor)):
        if nodes_length[i] > R_name[1]:
            nodes_r_st[i]= (intercept2 + slope2*nodes_length[i])
        else:
            nodes_r_st[i]= (intercept + slope*nodes_length[i])
    return {'connection':node_connection, 'organtype':nodes_organtype, 'resistance':nodes_r_st, 'unqiue':unq_cnt}


    #R_st is calculated based on the height (or the z axis value), 

    #calcuate the volumn of st


    # file-output.py


def write_PiafMunch_parameter(node_connection, nodes_organtype, nodes_r_st, unq_cnt): #function of creating parameters of PiafMunch
    Nt = len(nodes_organtype)-1 #total number of nodes, here the 0 and 1st node are the same, so minus 1
    Nc = len(node_connection) #total number of connections

    #condition = ==0
    #print(nodes_organtype[:,0]) #some nodes are not correct, here is the revise
    nodes_organtype[1,2] =2
    nodes_organtype[1,1] =1
    node_connection[0,0] =1
    #all the leaf source segment nodes
    N1L_node = nodes_organtype[(nodes_organtype[:,0] >2 ) & ((nodes_organtype[:,1] ==8)|(nodes_organtype[:,1] ==4)  ) & (nodes_organtype[:,2] ==1 )]
    #all the leaf source segments
    N1L_c_nd = list(range(1, len(N1L_node)+1))
    N1L_conn = list(range(1, len(N1L_node)+1))
    for i in range(len(N1L_node)):
        N1L_c_nd[i] = node_connection[(node_connection[:,0] == N1L_node[i][0])]
        N1L_conn[i] = np.where( node_connection[:,1] == N1L_c_nd[i][0,1])[-1]


    N1R_node = nodes_organtype[(nodes_organtype[:,0] >2 ) & (nodes_organtype[:,1] ==2 ) & (nodes_organtype[:,2] ==1 )]
    N1R_c_nd = list(range(1, len(N1R_node)+1))
    N1R_conn = list(range(1, len(N1R_node)+1))
    for i in range(len(N1R_node)):
        N1R_c_nd[i] = node_connection[(node_connection[:,1] == N1R_node[i][0])]
        N1R_conn[i] = np.where( node_connection[:,1] == N1R_c_nd[i][0,1])[0]
    N1R_r_abs = 1e-025

    N2_node = nodes_organtype[ (nodes_organtype[:,2] ==2 )]



    ################## Nodes With 2 Connections #########################
    N2_c_nd_1 = list(range(1, len(N2_node)+1))
    N2_conn_1 = list(range(1, len(N2_node)+1))
    for i in range(len(N2_node)):
        N2_c_nd_1[i] = node_connection[(node_connection[:,0] == N2_node[i,0])]
        N2_conn_1[i] = np.where( node_connection[:,0] == N2_c_nd_1[i][0,0])[0]


    #temp= N2_c_nd_1[0]
    #temp=temp[::-1] 
    #print(temp )

    N2_c_nd_2 = list(range(1, len(N2_node)+1))
    N2_conn_2 = list(range(1, len(N2_node)+1))

    for i in range(len(N2_node)):
        N2_c_nd_2[i] = node_connection[(node_connection[:,1] == N2_node[i,0])]
        N2_conn_2[i] = np.where( node_connection[:,1] == N2_c_nd_2[i][0][1])[0]

        

    ################## Nodes With 3 Connections #########################
    N3_node = nodes_organtype[ (nodes_organtype[:,2] ==3 )]
    #print(N3_node)

    N3_c_nd_1 = list(range(1, len(N3_node)+1))
    N3_conn_1 = list(range(1, len(N3_node)+1))
    for i in range(len(N3_node)):
        if N3_node[i,0] == N3_node[i-1,0]+1 : # some n3 nodes are inter connected
            N3_c_nd_1[i] = node_connection[(node_connection[:,1] == N3_node[i,0])][0][0]
            N3_conn_1[i] = 0-np.where( node_connection[:,0] ==  node_connection[(node_connection[:,0] == N3_node[i-1,0])][1][1])[0][0]
        elif  N3_node[i,1] == 2 or N3_node[i,1] == 1:
            N3_c_nd_1[i] = node_connection[(node_connection[:,1] == N3_node[i,0])][0][0]
            N3_conn_1[i] = 0-np.where( node_connection[:,0] == N3_c_nd_1[i])[0][0]-1
        
        else:
            N3_c_nd_1[i] = node_connection[(node_connection[:,0] == N3_node[i,0])][0][1]
            N3_conn_1[i] = np.where( node_connection[:,1] == N3_c_nd_1[i])[0][0]+1

    N3_c_nd_2 = list(range(1, len(N3_node)+1))
    N3_conn_2 = list(range(1, len(N3_node)+1))
    for i in range(len(N3_node)):
        if N3_node[i,1] == 2 or N3_node[i,1] == 1 and N3_node[i,0] != N3_node[i-1,0]+1:
            N3_c_nd_2[i] = node_connection[(node_connection[:,0] == N3_node[i,0])][0][1]
            N3_conn_2[i] = np.where( node_connection[:,1] == N3_c_nd_2[i])[0][0]+1
        else:
            N3_c_nd_2[i] = node_connection[(node_connection[:,1] == N3_node[i,0])][0][0]
            N3_conn_2[i] = 0-np.where( node_connection[:,0] == N3_c_nd_2[i])[0][0]-1

    N3_c_nd_3 = list(range(1, len(N3_node)+1))
    N3_conn_3 = list(range(1, len(N3_node)+1))
    for i in range(len(N3_node)):
        if N3_node[i,1] == 2 or N3_node[i,1] == 1 and N3_node[i,0] != N3_node[i-1,0]+1:
            N3_c_nd_3[i] = node_connection[(node_connection[:,0] == N3_node[i,0])][1][1]
            N3_conn_3[i] = np.where( node_connection[:,1] == N3_c_nd_3[i])[0][0]+1
        else:
            N3_c_nd_3[i] = node_connection[(node_connection[:,1] == N3_node[i,0])][1][0]
            N3_conn_3[i] = 0-np.where( node_connection[:,0] == N3_c_nd_3[i])[0][0]-1


    N4_node = nodes_organtype[ (nodes_organtype[:,2] ==4 )]
    #print(N4_node)

    N4_c_nd_1 = list(range(1, len(N4_node)+1))
    N4_conn_1 = list(range(1, len(N4_node)+1))
    for i in range(len(N4_node)):
        N4_c_nd_1[i] = node_connection[(node_connection[:,0] == N4_node[i,0])][0][1]
        N4_conn_1[i] = np.where( node_connection[:,1] == N4_c_nd_1[i])[0][0]+1

    #print(N4_c_nd_1)
    #print(N4_conn_1)
    #for i in range(len(N4_node)):
    #    print(N4_c_nd_1[i])
    #for i in range(len(N4_node)):
    #    print(N4_conn_1[i][0]+1)
        
        
    N4_c_nd_2 = list(range(1, len(N4_node)+1))
    N4_conn_2 = list(range(1, len(N4_node)+1))
    N4_conn_3 = list(range(1, len(N4_node)+1))
    N4_conn_4 = list(range(1, len(N4_node)+1))
    for i in range(len(N4_node)):
        N4_c_nd_2[i] = node_connection[(node_connection[:,1] == N4_node[i,0])]
        N4_conn_2[i] = np.where( node_connection[:,0] == N4_c_nd_2[i][0][0])[0]
        N4_conn_3[i] = np.where( node_connection[:,0] == N4_c_nd_2[i][1][0])[0]
        N4_conn_4[i] = np.where( node_connection[:,0] == N4_c_nd_2[i][2][0])[0]




    #'******** CARBON Lateral FLUX - RELATED PARAMETERS *********\n'
    #initialization of the parameters
    kML = np.zeros(len(nodes_organtype))
    vML = np.zeros(len(nodes_organtype))
    kMU = np.zeros(len(nodes_organtype))
    vMU = np.zeros(len(nodes_organtype))
    kMParMb = np.zeros(len(nodes_organtype))
    vMParMb = np.zeros(len(nodes_organtype))
    kM = np.zeros(len(nodes_organtype)) #kinetic parameter / Michaelis - starch Synthesis
    Vmax = np.zeros(len(nodes_organtype)) # kinetic parameter / starch Synthesis
    C_targ = np.zeros(len(nodes_organtype)) #kinetic parameter / starch/sugar equilibrium. (regul. par. sugar conc.)             (mmol / ml)
    kHyd = np.zeros(len(nodes_organtype))
    k1 = np.zeros(len(nodes_organtype))
    k2 = np.zeros(len(nodes_organtype))
    k3 = np.zeros(len(nodes_organtype))
    StructC = np.zeros(len(nodes_organtype))
    vol_ST = np.zeros(len(nodes_organtype))
    volPhlApo = np.zeros(len(nodes_organtype))
    volParApo = np.zeros(len(nodes_organtype))
    k_Lockhart = np.zeros(len(nodes_organtype))
    P_thr = np.zeros(len(nodes_organtype))
    vol_Sympl_max = np.zeros(len(nodes_organtype))

    r_Xyl = np.full(len(nodes_organtype), 0.0005)
    r_ST = np.full(len(nodes_organtype), nodes_r_st)       #automatically assign the sieve tube resistance calculated based on 
    r_Trsv = np.full(len(nodes_organtype), 100)
    r_PhlMb = np.full(len(nodes_organtype), 135.785)
    r_ParMb = np.full(len(nodes_organtype), 1e+025)
    r_Apo = np.full(len(nodes_organtype), 1e+025)
    r_Sympl = np.full(len(nodes_organtype), 1e+025)




    for i in range(len(nodes_organtype)): #given different value based on whether it is source, sink or connection
        if (nodes_organtype[i,1] == 8 or nodes_organtype[i,1] == 4) and nodes_organtype[i,2] == 1: #all the sources       
            kML[i]     =  1e-100
            vML[i]        = 0.000143136      #kinetic parameter / phloem loading (mmol /h) different in source, sink or connection of piafmunch2 oringinal value is 0.000143136 
            kMU[i]        = 10e-100     #different in source, sink or connection of piafmunch2
            vMU[i]        = 0      #different in source, sink or connection of piafmunch2
            kMParMb[i]    = 1
            vMParMb[i]    = 0
            kM[i]         = 1e-100
            Vmax[i]       = 0
            C_targ[i]     = 0.17      #different in source, sink or connection of piafmunch2
            kHyd[i]       = 0
            k1[i]         = 0
            k2[i]         = 0     #
            k3[i]         = 0
            StructC[i]    = 0      #different in source, sink or connection of piafmunch2
            vol_ST[i]     = 2.6e-05
            volPhlApo[i]  = 2.6e-05
            volParApo[i]  = 2.6e-05
            k_Lockhart[i] = 0
            P_thr[i]      = 1
            vol_Sympl_max[i] = 0.00018
        elif nodes_organtype[i,0]>0 and nodes_organtype[i,1]==2 and nodes_organtype[i,2]==1:   #all the sinks  
            kML[i]     =  1e-100
            vML[i]        = 0      #different in source, sink or connection of piafmunch2
            kMU[i]        =   1e+99      #different in source, sink or connection of piafmunch2 default 1e+99
            vMU[i]        = 2.82627e+95      #different in source, sink or connection of piafmunch2 default is 2.82627e+95
            kMParMb[i]    = 1
            vMParMb[i]    = 0
            kM[i]         = 1e-100
            Vmax[i]       = 0
            C_targ[i]     =     0.1     #different in source, sink or connection of piafmunch2
            kHyd[i]       = 0
            k1[i]         = 0
            k2[i]         = 0     #manually set it to 0.4
            k3[i]         = 0
            StructC[i]    = 1      #different in source, sink or connection of piafmunch2
            vol_ST[i]     = 2.6e-05
            volPhlApo[i]  = 2.6e-05
            volParApo[i]  = 2.6e-05
            k_Lockhart[i] = 0
            P_thr[i]      = 1
            vol_Sympl_max[i] = 0.00018        
        elif nodes_organtype[i,2]!=1: #all other connections other than source and sink
            kML[i]     =  1e-100
            vML[i]        = 0      #different in source, sink or connection of piafmunch2
            kMU[i]        = 1e-100      #different in source, sink or connection of piafmunch2
            vMU[i]        = 0      #different in source, sink or connection of piafmunch2
            kMParMb[i]    = 1
            vMParMb[i]    = 0
            kM[i]         = 1e-100
            Vmax[i]       = 0
            C_targ[i]     = 0.1      #different in source, sink or connection of piafmunch2
            kHyd[i]       = 0
            k1[i]         = 0
            k2[i]         = 0
            k3[i]         = 0
            StructC[i]    = 1      #different in source, sink or connection of piafmunch2
            vol_ST[i]     = 2.6e-05 #ml
            volPhlApo[i]  = 2.6e-05
            volParApo[i]  = 2.6e-05
            k_Lockhart[i] = 0
            P_thr[i]      = 1
            vol_Sympl_max[i] = 0.00018



    #'******** INITIAL VALUES *********\n'
    #initialization of the parameters
    Q_ST = np.full(len(nodes_organtype), 0)
    Q_Sympl = np.full(len(nodes_organtype), 4.4e-006)
    Starch = np.full(len(nodes_organtype), 1)
    Q_PhlApo = np.full(len(nodes_organtype), 4.4e-006)
    Q_ParApo = np.full(len(nodes_organtype), 4.4e-006)
    Tr_Q_ST = np.full(len(nodes_organtype), 0)
    Tr_Q_Sympl = np.full(len(nodes_organtype), 4.4e-006)
    Tr_Starch = np.full(len(nodes_organtype), 1)
    Tr_Q_PhlApo = np.full(len(nodes_organtype), 0)
    Tr_Q_ParApo = np.full(len(nodes_organtype), 0)
    vol_Sympl = np.full(len(nodes_organtype), 2.6e-005)

    #******** SIMULATION SOLVING PARAMETERS *********


    #'******** CARBON Lateral FLUX - RELATED PARAMETERS *********\n'
    #initialization of the parameters

    Q_ST_Abs = np.full(len(nodes_organtype), 2.6e-012)
    Q_Sympl_Abs = np.full(len(nodes_organtype), 1e-007)
    Starch_Abs = np.full(len(nodes_organtype), 1e-007)
    Q_PhlApo_Abs = np.full(len(nodes_organtype), 2.6e-012)
    Q_ParApo_Abs = np.full(len(nodes_organtype), 1e-007)
    Tr_Q_ST_Abs = np.full(len(nodes_organtype), 2.6e-012)
    Tr_Q_Sympl_Abs = np.full(len(nodes_organtype), 1e-007)
    Tr_Starch_Abs = np.full(len(nodes_organtype), 1e-007)
    Tr_Q_PhlApo_Abs = np.full(len(nodes_organtype), 2.6e-012)
    Tr_Q_ParApo_Abs = np.full(len(nodes_organtype), 1e-007)
    vol_Sympl_Abs = np.full(len(nodes_organtype), 2.6e-012)



    for i in range(len(nodes_organtype)): #given different value based on whether it is source, sink or connection
        if (nodes_organtype[i,1] == 8 or nodes_organtype[i,1] == 4) and nodes_organtype[i,2] == 1: #all the sources       
            Q_ST_Abs[i] =  1e-015
            Q_Sympl_Abs[i] =  1e-015
            Starch_Abs[i] =  1e-012
            Q_PhlApo_Abs[i] =  1e-015
            Q_ParApo_Abs[i] =  1e-015
            Tr_Q_ST_Abs[i] =  1e-012
            Tr_Q_Sympl_Abs[i] =  1e-012
            Tr_Starch_Abs[i] =  1e-012
            Tr_Q_PhlApo_Abs[i] =  1e-015
            Tr_Q_ParApo_Abs[i] =  1e-015
            vol_Sympl_Abs[i] = 1e-012


# Creat small functions for assign different speed or k

    def assign_source_loading_speed( nodes_organtype, value= 0.000143136):
        for i in range(len(nodes_organtype)): #given different value based on whether it is source, sink or connection
            if (nodes_organtype[i,1] == 8 or nodes_organtype[i,1] == 4) and nodes_organtype[i,2] == 1: #all the sources       
                vML[i]        = value      #kinetic parameter / phloem loading (mmol /h) different in source, sink or connection of piafmunch2 oringinal value is 0.000143136 
        return vML;

    def assign_source_loading_k( value= 1e-100):
        for i in range(len(nodes_organtype)): #given different value based on whether it is source, sink or connection
            if (nodes_organtype[i,1] == 8 or nodes_organtype[i,1] == 4) and nodes_organtype[i,2] == 1: #all the sources       
                kML[i]        = value      #kinetic parameter /  phloem loading (mmol /h) different in source, sink or connection of piafmunch2 oringinal value is 0.000143136 
        return kML;


    def assign_sink_unloading_speed(nodes_organtype,  value= 2.82627e+95):
        for i in range(len(nodes_organtype)): #given different value based on whether it is source, sink or connection
            if nodes_organtype[i,0]>0 and nodes_organtype[i,1]==2 and nodes_organtype[i,2]==1:   #all the sinks        
                vMU[i]        = value      #kinetic parameter / phloem loading (mmol /h) different in source, sink or connection of piafmunch2 oringinal value is 0.000143136 
        return vMU;

    def assign_sink_unloading_k( value= 1e+99):
        for i in range(len(nodes_organtype)): #given different value based on whether it is source, sink or connection
            if nodes_organtype[i,0]>0 and nodes_organtype[i,1]==2 and nodes_organtype[i,2]==1:   #all the sinks        
                kMU[i]        = value      #kinetic parameter / phloem loading (mmol /h) different in source, sink or connection of piafmunch2 oringinal value is 0.000143136 
        return kMU;

    def assign_resistance( value):
        r_ST = np.full(len(nodes_organtype), value) 
        return r_ST;




    def create_piafmunch_parameter( unq_cnt,  name= "mg_low1.ini" , end_time = "100"  ): #function of writing parameters of PiafMunch
        f = open(name,'w')
        f.write('******** DESCRIPTION OF ARCHITECTURE *********\n\n')

        f.write("Total number of Nodes : {0} = {1}\n".format('Nt', Nt))
        f.write("number of Internode Connections : {0} = {1}\n\n".format('Nc', Nc))

        f.write("Nodes Of Connectivity Order 1, Transpiring Leaf Ends : {0} = {1}\n".format('N1L', len(N1L_node)))
        f.write("{:s}  {:s}  {:s}\n".format('node#','c.node','conn.#'))
        for i in range(len(N1L_node)):
            f.write("{:d}  {:d}  {:d}\n\n".format(N1L_node[i][0],N1L_c_nd[i][0,1],(N1L_conn[i][-1]+1)))





        f.write("Nodes Of Connectivity Order 1, Absorbing Root Ends : {0} = {1}\n\n".format('N1R', len(N1R_node)))
        f.write("{:s}  {:s}  {:s}  {:s}\n".format('node#','c.node','conn.#','r_abs'))
        for i in range(len(N1R_node)):
            f.write("{:d}  {:d}  {:d} {:e}\n\n".format(N1R_node[i][0], N1R_c_nd[i][0,0],0- N1R_conn[i][0]-1, 1e-025))

        f.write('Nodes Of Connectivity Order 2 :  {0} = {1}\n\n' .format('N2', len(N2_node) ))
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','c.nd.1','conn.1','c.nd.2','conn.2'))
        for i in range(len(N2_node)):
            f.write("{:d}  {:d}  {:d} {:d} {:d} \n".format(N2_node[i][0],N2_c_nd_1[i][0,1],N2_conn_1[i][0]+1,N2_c_nd_2[i][0,0],0 -N2_conn_2[i][0]-1))    
        f.write('\n')

        f.write("Nodes Of Connectivity Order 3 :  {0} = {1}\n".format('N3', len(N3_node)))
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','c.nd.1','conn.1','c.nd.2','conn.2','c.nd.3','conn.3'))
        for i in range(len(N3_node)):
            f.write("{:d}  {:d}  {:d} {:d} {:d} {:d} {:d}\n".format(N3_node[i][0],N3_c_nd_1[i],N3_conn_1[i],N3_c_nd_2[i],N3_conn_2[i],N3_c_nd_3[i],N3_conn_3[i]))    
        f.write('\n')




        f.write("Nodes Of Connectivity Order 4 :  {0} = {1}\n\n".format('N4', np.count_nonzero(unq_cnt == 4)))
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','c.nd.1','conn.1','c.nd.2','conn.2','c.nd.3','conn.3','c.nd.4','conn.4'))

        for i in range(len(N4_node)):
            f.write("{:d}  {:d}  {:d} {:d} {:d} {:d} {:d} {:d} {:d}\n".format(N4_node[i][0],N4_c_nd_1[i],N4_conn_1[i],N4_c_nd_2[i][0][0],0-N4_conn_2[i][0]-1,N4_c_nd_2[i][1][0],0-N4_conn_3[i][0]-1,N4_c_nd_2[i][2][0],0-N4_conn_4[i][0]-1))    
        f.write('\n')



        f.write("Nodes Of Connectivity Order 5 :  {0} = {1}\n\n".format('N5', np.count_nonzero(unq_cnt == 5)))
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n\n".format('node#','c.nd.1','conn.1','c.nd.2','conn.2','c.nd.3','conn.3','c.nd.4','conn.4','c.nd.5','conn.5'))
        f.write("Nodes Of Connectivity Order 6 :  {0} = {1}\n".format('N6', np.count_nonzero(unq_cnt == 6)))
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n\n".format('node#','c.nd.1','conn.1','c.nd.2','conn.2','c.nd.3','conn.3','c.nd.4','conn.4','c.nd.5','conn.5','c.nd.6','conn.6'))
        f.write("Nodes Of Connectivity Order 7 :  {0} = {1}\n".format('N7', np.count_nonzero(unq_cnt == 7)))
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n\n".format('node#','c.nd.1','conn.1','c.nd.2','conn.2','c.nd.3','conn.3','c.nd.4','conn.4','c.nd.5','conn.5','c.nd.6','conn.6','c.nd.7','conn.7'))
        f.write("Nodes Of Connectivity Order 8 :  {0} = {1}\n".format('N8', np.count_nonzero(unq_cnt == 8)))
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n\n".format('node#','c.nd.1','conn.1','c.nd.2','conn.2','c.nd.3','conn.3','c.nd.4','conn.4','c.nd.5','conn.5','c.nd.6','conn.6','c.nd.7','conn.7','c.nd.8','conn.8'))

        f.write('******** WATER FLUX - RELATED PARAMETERS *********\n\n')
        f.write("{:s}  {:s}  {:s} \n".format('T\xb0 K','visco=f(C)','NonLin.Psi+NonZeroSugarVol.'))
        f.write("{:s}  {:s}  {:s} \n".format('293','true','true'))

        f.write('InterNode Connections -- Axial Resistances (MPa h / ml) : Nc= {}\n' .format(len(node_connection)))
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','upfl.node','dnfl.node','r_Xyl','r_ST' ))
        for i in range(len(node_connection)):
            f.write("{:d}  {:d}  {:d} {:e} {:e}\n".format(i+1,node_connection[i,0],node_connection[i,1],r_Xyl[i],r_ST[i]))
        f.write('\n')

        f.write('Individual Node : Lateral Resistances (MPa h / ml)\n')
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','r_Trsv','r_PhlMb','r_ParMb','r_Apo', 'r_Sympl' ))
        for i in range(len(nodes_organtype)-1):
            f.write("{:d}  {:e}  {:e} {:e} {:e} {:e}\n".format(nodes_organtype[i+1][0],r_Trsv[i],r_PhlMb[i],r_ParMb[i],r_Apo[i],r_Sympl[i]))
        f.write('\n')


        f.write('******** CARBON Lateral FLUX - RELATED PARAMETERS *********\n')    
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s}  {:s} \n".format('node#','kML(M)','vML(mmol/h)','kMU(M)','vMU(mmol/h)', 'kMParMb(M)','vMParMb(mmol/h)','kM(M)','Vmax(M/h)','C_targ(M)','kHyd(h-1)','k1(h-1)','k2','k3(h-1)','StructC','vol_ST(ml)','volPhlApo,ml','volParApo,ml','k_Lockhart','P_thr(MPa)','vol_Sympl_max,ml' ))    
        for i in range(len(nodes_organtype)-1):
            f.write("{:d}  {:e}  {:e} {:e} {:e} {:e}  {:e}  {:e} {:e} {:e} {:e}  {:e}  {:e} {:e} {:e} {:e}  {:e}  {:e} {:e} {:e} {:e}\n"
            .format(nodes_organtype[i+1][0],kML[i+1], vML[i+1], kMU[i+1], vMU[i+1], kMParMb[i+1], vMParMb[i+1], kM[i+1], Vmax[i+1], C_targ[i+1], kHyd[i+1], 
            k1[i+1], k2[i+1], k3[i+1], StructC[i+1], vol_ST[i+1], volPhlApo[i+1], volParApo[i+1], k_Lockhart[i+1], P_thr[i+1], vol_Sympl_max[i+1]))
        f.write('\n')



        f.write('******** INITIAL VALUES  *********\n')
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s} {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','Q.ST(mmol)','Q.Sympl(mmol)','Starch','Q.PhlApo(mmol)', 'Q.ParApo(mmol)',
                                                    'Tr.Q.ST(mmol)','Tr.Q.Sympl(mmol)','Tr.Starch','Tr.Q.PhlApo(mmol)', 'Tr.Q.ParApo(mmol)','vol_Sympl(ml)' ))
        for i in range(len(nodes_organtype)-1):
            f.write("{:d}  {:e}  {:e} {:e} {:e} {:e}  {:e}  {:e} {:e} {:e} {:e} {:e}\n"
                    .format(nodes_organtype[i+1][0], Q_ST[i+1],Q_Sympl[i+1],Starch[i+1],Q_PhlApo[i+1],Q_ParApo[i+1], 
                            Tr_Q_ST[i+1],Tr_Q_Sympl[i+1],Tr_Starch[i+1],Tr_Q_PhlApo[i+1],Tr_Q_ParApo[i+1],vol_Sympl[i+1]))
        f.write('\n')    

        f.write('******** SIMULATION SOLVING PARAMETERS *********\n')
        f.write('{:s}  {:s}  {:s}   {:s}  {:s}\n'.format('StartTime','EndTime','OutputStep', 'TracerHalfLife','Rel_Tol'))
        f.write('{:s}  {:s}  {:s}   {:s}  {:s}\n'.format('0', end_time, '0.166667', '0.33967', '1e-007'))
        f.write('\n')     

        f.write('***Abs_Tols for individual  nodes ***')
        f.write("{:s}  {:s}  {:s}  {:s}  {:s}  {:s} {:s}  {:s}  {:s}  {:s}  {:s}  {:s}\n".format('node#','Q.ST(mmol)','Q.Sympl(mmol)','Starch','Q.PhlApo(mmol)', 'Q.ParApo(mmol)',
                                                    'Tr.Q.ST(mmol)','Tr.Q.Sympl(mmol)','Tr.Starch','Tr.Q.PhlApo(mmol)', 'Tr.Q.ParApo(mmol)','vol_Sympl(ml)' ))
        for i in range(len(nodes_organtype)-1):
            f.write("{:d}  {:e}  {:e} {:e} {:e} {:e}  {:e}  {:e} {:e} {:e} {:e} {:e}\n"
                    .format(nodes_organtype[i+1][0], Q_ST_Abs[i+1],Q_Sympl_Abs[i+1],Starch_Abs[i+1],Q_PhlApo_Abs[i+1],Q_ParApo_Abs[i+1], 
                            Tr_Q_ST_Abs[i+1],Tr_Q_Sympl_Abs[i+1],Tr_Starch_Abs[i+1],Tr_Q_PhlApo_Abs[i+1],Tr_Q_ParApo_Abs[i+1],vol_Sympl_Abs[i+1]))
        f.write('\n')     

        f.write('******** OUTPUT SETTINGS : INDIVIDUAL NODE - LATERAL FLUXES-RELATED VARIABLES *********\n')    
        f.write('Nodes selected for plotting  : nsp = {:d}\n' .format((len(nodes_organtype))-1 ))    
        for i in range(len(nodes_organtype)-1):
            f.write('{:d}\n'.format(nodes_organtype[i+1][0]))

        f.write('individual-Node-related variables selected for plotting : nvp = 1\n')    
        f.write('C_ST (mmol / ml)\n')        
        #f.write('JS_PhlMb (mmol / h)\n')        
        #f.write('JW_Trsv (ml / h)\n')        


        f.write('Nodes selected for saving  : nss = {:d}\n' .format((len(nodes_organtype))-1 ))    
        for i in range(len(nodes_organtype)-1):
            f.write('{:d}\n'.format(nodes_organtype[i+1][0]))
        f.write('\n')    


        f.write('individual-Node-related variables selected for saving : nvs = 30\n')    
        f.write('''
        C_ApoUpflow (mmol / ml)
        C_ParApo (mmol / ml)
        C_PhlApo (mmol / ml)
        C_ST (mmol / ml)
        C_Sympl (mmol / ml)
        C_SymplUpflow (mmol / ml)
        JS_Apo (mmol / h)
        JS_ParMb (mmol / h)
        JS_PhlMb (mmol / h)
        JS_Sympl (mmol / h)
        JW_Apo (ml / h)
        JW_ParMb (ml / h)
        JW_Sympl (ml / h)
        JW_Trsv (ml / h)
        P_PhlApo (MPa)
        P_ST (MPa)
        P_ST_dot (MPa / h)
        P_Sympl (MPa)
        P_Sympl_dot (MPa / h)
        P_Xyl (MPa)
        PsiSoil (MPa)
        Psi_ParApo (MPa)
        Psi_PhlApo (MPa)
        Psi_ST (MPa)
        Q_PhlApo (mmol)
        Q_PhlApo_dot (mmol / h)
        Q_Sympl_dot (mmol / h)
        Transpirat (ml / h)
        vol_Sympl (ml)
        vol_Sympl_dot (ml / h)
        StarchSyn (mmol eq. Glu / h)
        Starch (mmol eq. Glu)
        ''')

        f.write('******** OUTPUT SETTINGS : INTERNODE CONNECTION - AXIAL FLUXES-RELATED VARIABLES *********\n')  
        f.write('node-to-node Fluxes selected for plotting  : fsp = {} \n' .format(len(node_connection)))
        for i in range(len(node_connection)):
            f.write("{:d}\n".format((i+1)))
        f.write('\n')

        f.write('node-to-node-Fluxes-related variables selected for plotting : fvp = 1\n')    
        #f.write('JS_ST (mmol / h)\n')        
        f.write('JW_ST (ml / h)\n')        


        f.write('node-to-node Fluxes selected for saving  : fss = {} \n' .format(len(node_connection)))
        for i in range(len(node_connection)):
            f.write("{:d}\n".format((i+1)))
        f.write('\n')

        f.write('node-to-node-Fluxes-related variables selected for saving : fvs = 4\n')    
        f.write('JS_ST (mmol / h)\n')        
        f.write('JW_ST (ml / h)\n')  
        f.write('C_Upflow (mmol / ml)\n')
        f.write('JW_Xyl (ml / h)\n')  

        print('output successful'.format())




        f.close()
        #return assign_source_loading_speed, assign_sink_unloading_speed, create_piafmunch_parameter;、
        #return create_piafmunch_parameter;
    return create_piafmunch_parameter;
