import vtk
import numpy as np
from rsml_writer import write_rsml as write_rsml2  # to write a rsml

""" 
VTK Tools, by Daniel Leitner (refurbished 6/2020) 

for vtk to numpy, and numpy to vtk conversions

reading: vtp, vtu 
writing: vtp, vtu, msh, dgf, rsml
"""


def vtk_points(p):
    """ creates vtkPoints from an numpy array """
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
    """ creates vtkCells from an numpy array """
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


def np_data(polydata, data_index = 0, cell = None):
    """ The cell or point data from vtkPolyData as numpy array 
    @param polydata    grid as a vtkPolyData object 
    @param data_index  index of the vtk cell or point data
    @param cell        True = cell data, False = point data, (default) None = auto detect
    @return cell or point data of data_index, a flag if it is cell data or not (interesting in case of auto detect)
    """
    if cell is None:  # auto detect
        try:
            data = polydata.GetCellData()
            p = data.GetArray(data_index)
            noa = p.GetNumberOfTuples()
            cell = True
        except:
            data = polydata.GetPointData()
            p = data.GetArray(data_index)
            noa = p.GetNumberOfTuples()
            cell = False
    else:  # use information
        if cell:
            data = polydata.GetCellData()
            p = data.GetArray(data_index)
            noa = p.GetNumberOfTuples()
        else:
            data = polydata.GetPointData()
            p = data.GetArray(data_index)
            noa = p.GetNumberOfTuples()
    p_ = np.ones(noa,)
    for i in range(0, noa):
        d = p.GetTuple(i)
        p_[i] = d[0]
    return p_, cell


def read_vtp(name):
    """ Opens a vtp and returns the vtkPolydata class, converts all cell data (additionally) to point data """
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
    """ Opens a vtu and returns the vtkUnstructuredGrid class """
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


def write_rsml(name, pd, meta, id_ind = 5):
    """ Writes a RMSL file from vtkPolyDat using rsml_reader.write_rsml """

    nodes = np_points(pd)
    try:
        segs = np_cells(pd)
    except:
        print("write_rsml error: write rsml expects a root system represented as line segments, propably the file consists of polylines")

    n = pd.GetPointData().GetNumberOfArrays()
    print("Node Data", n)
    node_data = np.zeros((n, nodes.shape[0]))
    for i in range(0, n):
        node_data [i, :], _ = np_data(pd, i, False)

    n = pd.GetCellData().GetNumberOfArrays()
    print("Cell Data", n)
    seg_data = np.zeros((n, segs.shape[0]))
    for i in range(0, n):
        print(i)
        seg_data[i, :], _ = np_data(pd, i, True)

    ids = np.array(seg_data[id_ind, :], dtype = int) + 2
    segs = np.array(segs, dtype = int)

    print("Orders:", ids)
    print("Segments")
    print(segs)
    print("Nodes")
    print(nodes)

    write_rsml2(name, [0], segs, ids, nodes, node_data, meta, Renumber = True)


def read3D_vtp_data(name, data_index = 0, cell = None):
    """ returns the cell or vertex data of vtp or vtu file and the corresponding coordinates
    @param name         filename of the vtp representing a root system
    @param cell         cell or point data 
    """
    if name.endswith('.vtp') or name.endswith('.VTP'):
        polydata = read_vtp(name)
    else:
        polydata = read_vtu(name)
    points = polydata.GetPoints()
    p_, cell = np_data(polydata, data_index, cell)
    if not cell:
        Np = polydata.GetNumberOfPoints()
        z_ = np.zeros((Np, 3))
        for i in range(0, Np):
            p = np.zeros(3,)
            points.GetPoint(i, p)
            z_[i, :] = p
        return p_, z_
    if cell:
        Nc = polydata.GetNumberOfLines()
        z_ = np.zeros((Nc, 3))
        for i in range(0, Nc):
            c = polydata.GetCell(i)
            ids = c.GetPointIds()
            p1 = np.zeros(3,)
            points.GetPoint(ids.GetId(0), p1)
            p2 = np.zeros(3,)
            points.GetPoint(ids.GetId(1), p2)
            z_[i, :] = 0.5 * (p1 + p2)
        return p_, z_


def read3D_vtp_data_parallel(prename, postname, n, data_index = 0):
    """ returns the cell or vertex data of parallel vtp files """
    z_ = np.ones((0, 3))
    d_ = np.ones(0,)
    for i in range(0, n):
        n_ = prename + "{:04d}-".format(i) + postname + ".vtu"
        print("opening: ", n_)
        d, z = read3D_vtp_data(n_, data_index, None)
        z_ = np.hstack((z_, z))
        d_ = np.hstack((d_, d))
    return d_, z_


def read3D_data(name, np = 1, data_index = 0):
    """ Opens a vtu (parallel or not) and returns the vtkUnstructuredGrid class """
    if np == 1:
        print("read3D_data: open single")
        return read3D_vtp_data(name + ".vtu", data_index, None)
    else:
        print("read3D_data: open multiple")
        return read3D_vtp_data_parallel("s{:04d}-p".format(np), name, np, data_index)

