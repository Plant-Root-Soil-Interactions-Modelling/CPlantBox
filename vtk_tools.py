import vtk
import numpy as np

# from rsml.rsml_writer import write_rsml as write_rsml2  # to write a rsml (no other dependencies)
# from rsml.rsml_writer import Metadata

""" 
VTK Tools, by Daniel Leitner (refurbished 6/2020) 

Used by CPlantBox and dumux-rosi 

for post-processing of vtp and vtu files
for simpler vtk to numpy, and numpy to vtk conversions

reading: vtp, vtu (see read3D_data)
writing: vtp, vtu, msh, dgf, rsml
"""


def vtk_points(p):
    """ creates vtkPoints from an numpy array @param p"""
    points = vtk.vtkPoints()
    points.SetData(vtk_data3(p))
    return points


def vtk_cells(t):
    """ creates vtkCells from an numpy array @param t, 
    currently vtkLine, vtkTriangle, andn vtkTetra are implemented"""
    cellArray = vtk.vtkCellArray()
    for vert in t:
        if t.shape[1] == 2:
            tetra = vtk.vtkLine()
        if t.shape[1] == 3:
            tetra = vtk.vtkTriangle()
        elif t.shape[1] == 4:
            tetra = vtk.vtkTetra()
        for i, v in enumerate(vert):
            tetra.GetPointIds().SetId(i, int(v))
        cellArray.InsertNextCell(tetra)
    return cellArray


def vtk_data(d):
    """ creates a vtkDataArray from an numpy array @param d, 
     
    see also: grid.GetCellData().SetScalars(vtk_data(celldata)), grid.GetCellData().AddArray(...)
    
    TODO vtkAbstractArray.SetComponentName
    TODO auto detect number of tuples
    """
    da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
    da.SetNumberOfComponents(1)  # number of components
    da.SetNumberOfTuples(len(d))
    for i in range(0, len(d)):
        da.InsertTuple1(i, d[i])
#         elif noc == 2: # number of components
#             da.InsertTuple2(i, d[i, 0], d[i, 1])
#         elif noc == 3:
#             da.InsertTuple3(i, d[i, 0], d[i, 1], d[i, 2])
#         elif noc == 4:
#             da.InsertTuple4(i, d[i, 0], d[i, 1], d[i, 2], d[i, 3])
    return da


def vtk_data3(p):
    """ creates vtkPoints from an numpy array @param p"""
    da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
    da.SetNumberOfComponents(3)  # vtk point dimension is always 3
    da.SetNumberOfTuples(p.shape[0])
    for i in range(0, p.shape[0]):
        if p.shape[1] == 2:
            da.InsertTuple3(i, p[i, 0], p[i, 1], 0.)
        elif p.shape[1] == 3:
            da.InsertTuple3(i, p[i, 0], p[i, 1], p[i, 2])
    return da


def np_convert(x):
    """ converts a list of something (e.g. list of 3 floats) to a numpy array (Nx3)"""
    return np.array(list(map(np.array, x)))  # is there a better way?


def np_points(pd):
    """ returns the points of vtkPolyData as numpy array (Nx3)"""
    Np = pd.GetNumberOfPoints()
    z_ = np.zeros((Np, 3))
    points = pd.GetPoints()
    for i in range(0, Np):
        p = np.zeros(3,)
        points.GetPoint(i, p)
        z_[i,:] = p
    return z_


def np_cells(pd):
    """ returns the cells of vtkPolyData as numpy array (NxNc) 
    all cells must have the same number of points, 
    therefore it will not work for polylines """
    Nc = pd.GetNumberOfCells()
    d = pd.GetCell(0).GetPointIds().GetNumberOfIds()
    z_ = np.zeros((Nc, d), dtype = np.int64)
    for i in range(0, Nc):
        p = np.zeros(d,)
        ids = pd.GetCell(i).GetPointIds()
        for j in range(0, d):
            p[j] = ids.GetId(j)
        z_[i,:] = p
    return z_


def np_data(pd, data_index = 0, cell = None):
    """ returns cell or point data from vtkPolyData as numpy array 
    @param pd          rid as a vtkPolyData object 
    @param data_index  index of the vtk cell or point data (default = 0)
    @param cell        True = cell data, False = point data, (default = None, is auto detect)
    
    @return p_, cell   cell or point data of data_index, 
                       a flag if it is cell data or not (interesting in case of auto detect)
    """
    if cell is None:  # auto detect
        try:
            data = pd.GetCellData()
            p = data.GetArray(data_index)
            noa = p.GetNumberOfTuples()
            cell = True
        except:
            data = pd.GetPointData()
            p = data.GetArray(data_index)
            noa = p.GetNumberOfTuples()
            cell = False
    else:  # use information
        if cell:
            data = pd.GetCellData()
            p = data.GetArray(data_index)
            noa = p.GetNumberOfTuples()
        else:
            data = pd.GetPointData()
            p = data.GetArray(data_index)
            noa = p.GetNumberOfTuples()
    p_ = np.ones(noa,)
    for i in range(0, noa):
        d = p.GetTuple(i)
        p_[i] = d[0]
    return p_, cell


def rebuild_grid(p, t):
    """ Deletes unused points and updates elements
    @param p         nodes Nx3
    @param t         cells (as node indices) NxNc
    """
    pp = np.zeros(p.shape[0], dtype = "bool")  # initially all are unused
    for t_ in t:  # find unused points
        for n in t_:
            pp[n] = 1  # used
    upi = np.argwhere(pp == 0)  # unused point indices
    for k in upi[::-1]:  # reverse
        for i, t_ in enumerate(t):  # update triangle indices
            for j, n in enumerate(t_):
                if n > k:
                    t[i][j] -= 1
    p = np.delete(p, upi, axis = 0)  # delete unused points
    return p, t


def snap_to_box(p, box, eps = 1e-6):
    """ Snap points to a bounding box 
    @param p          nodes Nx3
    @param box        bounding box [minx,miny,minz,maxx,maxy,maxz]
    @param eps        snap tolerance
    """
    for i, p_ in enumerate(p):
        for j in range(0, 3):
            if p_[j] < box[j] + eps:
                p[i, j] = box[j]
        for j in range(3, 6):
            if p_[j - 3] > box[j] - eps:
                p[i, j - 3] = box[j]
    return p


def read_vtp(name):
    """ opens a vtp and returns the vtkPolydata class, 
    additionally, converts all cell datato point data (to plot the root system radius must be point data)
    @param name         including the file extension ("name.vtp") 
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(name)
    reader.Update()
    pd = reader.GetOutput()
    c2p = vtk.vtkCellDataToPointData()  # set cell and point data
    c2p.SetPassCellData(True)
    c2p.SetInputData(pd)
    c2p.Update()
    return c2p.GetPolyDataOutput()


def read_vtu(name):
    """ opens a vtu and returns the vtkUnstructuredGrid class
    @param name         file name without file extension ("name")   
    """
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(name)
    reader.Update()
    ug = reader.GetOutput()
    return ug


def read_rect_vtu(name):
    """ opens a vtu and returns the vtkXMLImageDataReader class 
    @param name         file name without file extension ("name")       
    """
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(name)
    reader.Update()
    ug = reader.GetOutput()
    return ug


def write_msh(name, pd):
    """ writes a tetraedral .msh file including cell data from vtkPolyData 
    @param name          file name including the file extension ("name.msh")   
    @param pd            grid represented as vtkPolyData object
    """
    with open(name, "w") as f:
        # Init
        f.write("$MeshFormat\n")
        f.write("2.2 0 8\n")  # version, file-type (ascii=0), data-size
        f.write("$EndMeshFormat\n")
        # Nodes
        np_ = pd.GetNumberOfPoints()
        f.write("$Nodes\n")
        f.write("{:d}\n".format(np_))  # number of nodes
        for i in range(0, np_):
            p = pd.GetPoint(i)
            f.write('{:d} {:08.6f} {:08.6f} {:08.6f}\n'.format(i + 1, p[0], p[1], p[2]))  # node number, x, y, z
        f.write("$EndNodes\n")
        # Cells
        ind = np.zeros(4, dtype = int)
        nc = pd.GetNumberOfCells()
        cdata = pd.GetCellData()
        dn = cdata.GetNumberOfArrays()
        f.write("$Elements\n")
        f.write("{:d}\n".format(nc + 8))  # number of cells
        for i in range(0, 8):
            f.write("{:d} 15 1 1 {:d}\n".format(i + 1, i + 1))
        for i in range(0, nc):
            tetra = pd.GetCell(i)
            c = tetra.GetPointIds()
            f.write("{:d} 4 {:d} ".format(i + 1 + 8, dn))  # id, 4 = tetra
            if dn > 0:
                for j in range(0, dn):
                    cdataj = cdata.GetArray(j)
                    d = cdataj.GetTuple(i)
                    f.write("{:g} ".format(d[0]))
            for j in range(0, 4):
                ind[j] = c.GetId(j) + 1

            f.write("{:d} {:d} {:d} {:d}\n".format(ind[0], ind[1], ind[2], ind[3]))
        f.write("$EndElements\n")


def write_dgf(name, pd):
    """ writes a DGF file including cell and point data from vtkPolyData
    @param name          file name including the file extension ("name.dgf")   
    @param pd            grid represented as vtkPolyData object
    """
    file = open(name, "w")  # write dgf
    file.write("DGF\n")
    # vertex
    file.write('Vertex\n')
    Np = pd.GetNumberOfPoints()
    points = pd.GetPoints()
    pdata = pd.GetPointData()
    Npd = pdata.GetNumberOfArrays()
    file.write('parameters {:g}\n'.format(Npd))
    for i in range(0, Np):
        p = np.zeros(3,)
        points.GetPoint(i, p)
        file.write('{:g} {:g} {:g} '.format(p[0], p[1], p[2]))
        for j in range(0, Npd):
            pdataj = pdata.GetArray(j)
            d = pdataj.GetTuple(i)
            file.write('{:g} '.format(d[0]))
        file.write('\n')
    file.write('#\n');
    # cell
    file.write('Simplex\n');
    Nc = pd.GetNumberOfCells()
    cdata = pd.GetCellData()
    Ncd = cdata.GetNumberOfArrays()
    file.write('parameters {:g}\n'.format(Ncd))
    for i in range(0, Nc - 1):
        cpi = vtk.vtkIdList()
        pd.GetCellPoints(i, cpi)
        for j in range(0, cpi.GetNumberOfIds()):  # write cell ids
            file.write('{:g} '.format(cpi.GetId(j)))
        for j in range(0, Ncd):
            cdataj = cdata.GetArray(j)
            d = cdataj.GetTuple(i)
            file.write('{:g} '.format(d[0]))
        file.write('\n')
    # i dont know how to use the following in dumux
    file.write('#\n')
    file.write('BOUNDARYSEGMENTS\n')  # how do i get the boundary segments into DUMUX ?
    file.write('2 0\n')
    file.write('3 {:g}\n'.format(Np - 1))  # vertex id, but index starts with 0
    file.write('#\n');
    file.write('BOUNDARYDOMAIN\n')
    file.write('default 1\n');
    file.write('#\n')
    file.close()


def write_vtp(name, pd):
    """ Writes a VTP file including cell data from vtkPolyData 
    @param name          file name including the file extension ("name.vtp")   
    @param pd            grid represented as vtkPolyData object
    """
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(name)
    writer.SetInputData(pd)
    writer.Write()


def write_vtu(name, pd):
    """ Writes a VTU file
    @param name          file name including the file extension ("name.vtp")   
    @param pd            grid represented as vtkPolyData object
    """
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(name)
    writer.SetInputData(pd)
    writer.Write()


def write_rsml(filename, pd, id_ind, meta = None, axes = [0]):
    """ Writes a RMSL file from vtkPolyData using rsml_reader.write_rsml 
    uses RSML function tags to store all vtk node or cell data    
    @param filename      output file name
    @param pd            vtk poly data file consisting of line segments
    @param id_ind        data array of order, or branch number, or something const along each branch (>=-1)
    @param meta          MetaData object (optional), see rsml_writer.py
    @param axes          list of indices of the starting nodes
    """
    nodes = np_points(pd)
    try:
        segs = np_cells(pd)
    except:
        print("write_rsml error: write rsml expects a root system represented as line segments, propably the file consists of polylines")

    n0 = pd.GetCellData().GetNumberOfArrays()
    n1 = pd.GetPointData().GetNumberOfArrays()
    data = -1 * np.ones((n0, nodes.shape[0]))  # -1 indicates not set
    names = []

    vtk_cell_data = pd.GetCellData()
    print("Cell Data (translated to rsml functions)", n0)
    for i in range(0, n0):
        name = vtk_cell_data.GetArray(i).GetName()
        names.append(name)
        seg_data, _ = np_data(pd, i, True)
        for j in range(0, segs.shape[0]):  # convert point to cell data
            data[i, segs[j, 1]] = seg_data[j]
        if i == id_ind:
            ids = np.array(seg_data, dtype = int) + 2  # needs to be >0 for reconstruction!
        print("\t", name, "\t[", np.min(data[i, 1:]), ", ", np.max(data[i, 1:]), "]")

    for i in range(0, n0):  # the first nodes of base polylines need to be set
        for j in range(0, nodes.shape[0] - 1):
            if data[i, j] < 0:
                data[i, j] = data[i, j + 1]

    print("Node data (rsml functions)", n1)
    vtk_node_data = pd.GetPointData()
    for i in range(0, n1):
        name = vtk_node_data.GetArray(i).GetName()
        if name and (not name in names):  # only, if not available as cell data
            names.append(name)
            f, _ = np_data(pd, i, False)
            data = np.vstack(data, f)
            print("\t", name, "\t[", np.min(data[i,:]), ", ", np.max(data[i,:]), "]")

    print("\nReconstruct from:", names[id_ind])  # should be orders
    print("\nSegments", segs.shape)
    print("Nodes", nodes.shape)

    if not meta:
        meta = Metadata()
    meta.set_fun_names(names)

    write_rsml2(filename, axes, segs, ids, nodes, data, meta, Renumber = True)


def read3D_vtp_data(name, data_index = 0, cell = None):
    """ returns the cell or vertex data of vtp or vtu file and the corresponding coordinates
    @param name         filename of the vtp ("name.vtp") or vtu ("name" without file extension) 
    @param data_index   index of cell or point data 
    @param cell         cell data (True) or point data (False), or default auto detect (None)
    
    @return p_, z_      cell or point data p, at coordinates z (cell centers or grid coordinates)
    """
    if name.endswith('.vtp') or name.endswith('.VTP'):
        pd = read_vtp(name)
    else:
        pd = read_vtu(name)
    points = pd.GetPoints()
    p_, cell = np_data(pd, data_index, cell)
    if not cell:
        Np = pd.GetNumberOfPoints()
        z_ = np.zeros((Np, 3))
        for i in range(0, Np):
            p = np.zeros(3,)
            points.GetPoint(i, p)
            z_[i,:] = p
        return p_, z_
    if cell:
        Nc = pd.GetNumberOfCells()
        z_ = np.zeros((Nc, 3))
        for i in range(0, Nc):
            c = pd.GetCell(i)
            ids = c.GetPointIds()
            n = ids.GetNumberOfIds ()
            midz = 0.
            p1 = np.zeros(3,)
            for j in range(0, n):
                points.GetPoint(ids.GetId(j), p1)
                midz += p1 / n
            z_[i,:] = midz
        return p_, z_


def read3D_vtp_data_parallel(prename, postname, n, data_index = 0, cell = None):
    """ cell or vertex data of parallel DuMux vtp, or vtu files, called by read3D_data 
    @param prename         the first part of the filename 
    @param postname        the name of the DuMux simulation
    @param n               number of processes used for parallel computation
    @param data_index      index of cell or point data
    @param cell         cell data (True) or point data (False), or default auto detect (None)
    
    @return d_, z_         cell or point data p, at coordinates z (cell centers or grid coordinates)
    """
    z_ = np.ones((0, 3))
    d_ = np.ones(0,)
    for i in range(0, n):
        n_ = prename + "{:04d}-".format(i) + postname + ".vtu"
        print("opening: ", n_)
        d, z = read3D_vtp_data(n_, data_index, cell)
        z_ = np.vstack((z_, z))
        d_ = np.hstack((d_, d))
    return d_, z_


def read3D_data(name, np = 1, data_index = 0, cell = None):
    """ opens a vtp or vtu (parallel or not) 
        @param np            number of processes used in parallel computation
        @param data_index    index of cell or point data 
        @param cell         cell data (True) or point data (False), or default auto detect (None)
        
        @return p, z         cell or point data p, at coordinates z (cell centers or grid coordinates)    
    """
    if np == 1:
        print("read3D_data: open single")
        if name.endswith('.vtp') or name.endswith('.VTP'):
            return read3D_vtp_data(name, data_index, cell)
        else:
            return read3D_vtp_data(name + ".vtu", data_index, cell)
    else:
        print("read3D_data: open multiple")
        return read3D_vtp_data_parallel("s{:04d}-p".format(np), name, np, data_index)

