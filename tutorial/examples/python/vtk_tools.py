import vtk
import numpy as np
# from rsml_writer import write_rsml # to write a rsml

""" 
VTK Tools, by Daniel Leitner (refurbished 6/2020) 

for vtk to numpy, and numpy to vtk conversions

reading: vtp, vtu 
writing: vtp, vtu, msh, dgf, rsml
"""


def vtk_points(p):
    """ Creates vtkPoints from an numpy array """
    da = vtk.vtkDataArray.CreateDataArray(vtk.VTK_DOUBLE)
    da.SetNumberOfComponents(3)  # vtk point dimension is always 3
    da.SetNumberOfTuples(p.shape[0])
    for i in range(0, p.shape[0]):
        if p.shape[1] == 2:
            da.InsertTuple3(i, p[i, 0], p[i, 1], 0.)
        elif p.shape[1] == 3:
            da.InsertTuple3(i, p[i, 0], p[i, 1], p[i, 2])
    points = vtk.vtkPoints()
    points.SetData(da)
    return points


def vtk_cells(t):
    """ Creates vtkCells from an numpy array """
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
    """ Creates a vtkDataArray from an numpy array, usage 
    e.g. grid.GetCellData().SetScalars(vtk_data(celldata)), grid.GetCellData().AddArray(...)
    TODO vtkAbstractArray.SetComponentName
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


def np_convert(x):
    """ Converts a list of something (list of 3 floats) to a numpy array (Nx3)"""
    return np.array(list(map(np.array, x)))  # is there a better way?


def np_points(polydata):
    """ The points of vtkPolyData as numpy array """
    Np = polydata.GetNumberOfPoints()
    z_ = np.zeros((Np, 3))
    points = polydata.GetPoints()
    for i in range(0, Np):
        p = np.zeros(3,)
        points.GetPoint(i, p)
        z_[i, :] = p
    return z_


def np_cells(polydata):
    """ The cells of vtkPolyData as numpy array  """
    Nc = polydata.GetNumberOfCells()
    d = polydata.GetCell(0).GetPointIds().GetNumberOfIds()
    z_ = np.zeros((Nc, d))
    for i in range(0, Nc):
        p = np.zeros(d,)
        ids = polydata.GetCell(i).GetPointIds()
        for j in range(0, d):
            p[j] = ids.GetId(j)
        z_[i, :] = p
    return z_


def np_data(polydata, data_index = 0, cell = True):
    """ The cell or vertex data from vtkPolyData as numpy array """
    if cell:
        data = polydata.GetCellData()
    else:
        data = polydata.GetPointData()
    p = data.GetArray(data_index)
    noa = p.GetNumberOfTuples()
    p_ = np.ones(noa,)
    for i in range(0, noa):
        d = p.GetTuple(i)
        p_[i] = d[0]
    return p_


def rebuild_grid(p, t):
    """ Deletes unused points and updates elements """
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
    """ Snap points to a bounding box """
    for i, p_ in enumerate(p):
        for j in range(0, 3):
            if p_[j] < box[j] + eps:
                p[i, j] = box[j]
        for j in range(3, 6):
            if p_[j - 3] > box[j] - eps:
                p[i, j - 3] = box[j]
    return p


def read_vtp(name):
    """ Opens a vtp and returns the vtkPolydata class """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(name)
    reader.Update()
    polydata = reader.GetOutput()
    return polydata


def read_vtu(name):
    """ Opens a vtp and returns the vtkUnstructuredGrid class """
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(name)
    reader.Update()
    ug = reader.GetOutput()
    return ug


def read_rect_vtu(name):
    """ Opens a vtp and returns the vtkUnstructuredGrid class """
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(name)
    reader.Update()
    ug = reader.GetOutput()
    return ug


def write_msh(name, pd):
    """ Writes a tetraedral .msh file including cell data from vtkPolyData """
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
    """ Writes a DGF file including cell and point data from vtkPolyData """
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
    """ Writes a VTP file including cell data from vtkPolyData """
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(name)
    writer.SetInputData(pd)
    writer.Write()


def write_vtu(name, grid):
    """ Writes a VTU file """
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(name)
    writer.SetInputData(grid)
    writer.Write()


def write_rsml(name, pd, meta):
    """ Writes a RMSL file from vtkPolyDat using rsml_reader.write_rsml """

    nodes = np_points(pd)
    segs = np_cells(pd)

    n = pd.GetPointData().GetNumberOfArrays()
    print("Node Data", n)
    node_data = np.zeros((n, nodes.shape[0]))
    for i in range(0, n):
        node_data [i, :] = np_data(pd, i, False)
    n = pd.GetCellData().GetNumberOfArrays()

    print("Cell Data", n)
    seg_data = np.zeros((n, segs.shape[0]))
    for i in range(0, n):
        seg_data[i, :] = np_data(pd, i, True)

    types = np.array(seg_data[3, :], dtype = int) + 2
    segs = np.array(segs, dtype = int)

    print("Orders:", types)
    print("Segments")
    print(segs)
    print("Nodes")
    print(nodes)

    write_rsml(name, [0], segs, types  , nodes, node_data, meta, Renumber = False)

