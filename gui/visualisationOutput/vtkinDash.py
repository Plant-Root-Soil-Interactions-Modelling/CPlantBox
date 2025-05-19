import os
import dash
from dash import html
import numpy as np
import dash_vtk
from dash_vtk.utils import to_volume_state, to_mesh_state
from vtk.util.numpy_support import vtk_to_numpy
import vtk


def uniform_grid(min_, max_, res):
    """ Creates an uniform grid
    @param min_    minimum of bounding rectangle
    @param max_    maximum of bounding rectangle
    @param res_    cell resolution
    @return A vtkUniformGrid
    """
    grid = vtk.vtkUniformGrid()
    grid.SetDimensions(int(res[0]) + 1, int(res[1]) + 1, int(res[2]) + 1)  # cells to corner points
    grid.SetOrigin(min_[0], min_[1], min_[2])
    s = (max_ - min_) / res
    grid.SetSpacing(s[0], s[1], s[2])
    return grid
def get_bounds_from_volume_state(volume_state):
    image = volume_state["image"]
    origin = np.array(image["origin"])
    spacing = np.array(image["spacing"])
    dimensions = np.array(image["dimensions"])

    # The bounds go from origin to origin + spacing * (dim - 1)
    max_bound = origin + spacing * (dimensions - 1)
    min_bound = origin  # just the origin

    return min_bound, max_bound
    
def plot_mesh_cuts(grid, array_name, nz = 3):
    """ plots orthogonal nz vertical cuts z[:-1] (xy-planes), with z = linspace(min_z, max_z, nz+1),
    and two additonal sclices at x=0 (yz-plane), y=0 (xz-plane)
    @param grid         some vtk grid (structured or unstructured)
    @param array_name       parameter to visualize
    @param nz           number of vertical slices
    """

    eps = 1.e-2
    planes = []  # create the cut planes
    bounds = grid.GetBounds()
    z = np.linspace(bounds[4] + eps, bounds[5], nz + 1)
    for i in range(0, nz):  # z-slices (implicit functions)
        p = vtk.vtkPlane()
        p.SetOrigin(0, 0, z[i])
        p.SetNormal(0, 0, 1)
        planes.append(p)
    for n in [(1, 0, 0), (0, 1, 0)]:
        p = vtk.vtkPlane()
        p.SetOrigin(bounds[0] + eps, bounds[2] + eps, bounds[4])
        p.SetNormal(n[0], n[1], n[2])
        planes.append(p)
    output_states = []
    for p in planes:
        cutter = vtk.vtkCutter()
        cutter.SetInputData(grid)

        cutter.SetCutFunction(p)
        cutter.Update()
        cut_output = cutter.GetOutput()
        # Ensure scalar array is active
        cut_output.GetCellData().SetActiveScalars(array_name)

        mesh_state = to_mesh_state(cut_output, field_to_keep=array_name)
        output_states.append(mesh_state)

    return output_states
    
# Data file path
demo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
head_vti = os.path.join(demo_dir,"visualisationOutput", "data", "head.vti")


# Load dataset from dist
reader = vtk.vtkXMLImageDataReader()
reader.SetFileName(head_vti)
reader.Update()
image_data = reader.GetOutput()
cell_data = image_data.GetCellData()
pname = "pressure head"
array = cell_data.GetArray(pname)
data_array = vtk_to_numpy(array).flatten()
data_range = [float(data_array.min()), float(data_array.max())]

if False:
    volume_state = to_volume_state(image_data)

    dimensions = volume_state["image"]["dimensions"]
    cell_count = [dimensions[0] - 1 , dimensions[1] - 1 , dimensions[2] - 1]
    min_bound, max_bound = get_bounds_from_volume_state(volume_state)
    ugrid = uniform_grid(min_bound, max_bound , cell_count)


    ugrid.GetCellData().AddArray(array)
    ugrid.GetCellData().SetActiveScalars(pname)
    volume_state = to_volume_state(ugrid)
    # Add the scalar field manually


    volume_state["field"] = {
        "name": pname,
        "location": "cellData",
        "data": data_array.tolist(),
        "dataRange": data_range,
        "lookupTable": {
            "mode": "opacity",
            "preset": "erdc_rainbow_bright",  # still required for colors
            "range": data_range,
            "mappingRange": data_range,
            "opacityPoints": [
                [data_range[0], 0.0],
                [data_range[1], 0.5]
            ]
        }
    }

if False:
    vtk_view = dash_vtk.View(
            children=[
                dash_vtk.VolumeRepresentation(
                    children=[
                        dash_vtk.Volume(state=volume_state),
                    ],
                    colorMapPreset="erdc_rainbow_bright",
                    colorDataRange=[
                        float(array.GetRange()[0]),
                        float(array.GetRange()[1])
                    ],
                )
            ]
    )
else:
    cut_states =plot_mesh_cuts(image_data, pname, nz = 3)
    vtk_view = dash_vtk.View(
    children=[
        dash_vtk.GeometryRepresentation(
            children=[
                dash_vtk.Mesh(state=state)
            ],
            colorMapPreset="erdc_rainbow_bright",
            colorDataRange=data_range,
        ) for state in cut_states
    ]
    )


app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    style={"height": "calc(100vh - 16px)", "width": "100%"},
    children=[html.Div(vtk_view, style={"height": "100%", "width": "100%"})],
)

if __name__ == "__main__":
    app.run(debug=True, port=8051)
