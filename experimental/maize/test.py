import numpy as np
import vtk


def create_maize_kernel_mesh(
    length=14.26,
    w_bottom=9.61,
    w_top=4.57,
    thickness=4.53,
    l_pedicel=1.8,
    n_x=120,
    n_theta=80,
    exponent=3.5,
):
    """Generates a 3D VTK PolyData surface mesh representing the single maize kernel model (Fig. 5b).

    Dimensions (in mm):
      - length: 14.26 mm (longitudinal length of main body)
      - w_bottom: 9.61 mm (crown width)
      - w_top: 4.57 mm (width at pedicel base)
      - thickness: 4.53 mm (ventral-dorsal thickness)
      - l_pedicel: 1.8 mm (pedicel tip protrusion length)
    """
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()

    theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)

    # -------------------------------------------------------------
    # 1. Main Kernel Body Mesh (x from 0 to length)
    # -------------------------------------------------------------
    n_body = int(n_x * 0.85)
    x_body = np.linspace(0, length, n_body)
    grid_pts = []

    for x in x_body:
        zeta = x / length  # Normalized position [0, 1]

        if zeta < 0.86:
            # Tapering from top width (4.57 mm) to crown width (9.61 mm)
            w_x = w_top + (w_bottom - w_top) * (zeta / 0.86) ** 0.92
            t_x = thickness * (1.0 - 0.12 * (zeta - 0.45) ** 2)
        else:
            # Crown dome at top end (with subtle dent curvature)
            cap_ratio = (zeta - 0.86) / 0.14
            cap_factor = np.sqrt(max(0.0, 1.0 - cap_ratio**2))
            w_x = w_bottom * cap_factor
            t_x = thickness * 0.98 * cap_factor

        a_y = max(w_x / 2.0, 0.01)
        b_z = max(t_x / 2.0, 0.01)

        row = []
        for th in theta:
            cos_t = np.cos(th)
            sin_t = np.sin(th)

            # Superelliptic cross-section for smooth rounded-trapezoidal shape
            y = a_y * np.sign(cos_t) * (np.abs(cos_t) ** (2.0 / exponent))
            z = b_z * np.sign(sin_t) * (np.abs(sin_t) ** (2.0 / exponent))

            # Crown center indentation
            if zeta > 0.88 and np.abs(cos_t) < 0.3:
                z *= 1.0 - 0.08 * np.sin(np.pi * (zeta - 0.88) / 0.12)

            pt_id = points.InsertNextPoint(x, y, z)
            row.append(pt_id)
        grid_pts.append(row)

    # -------------------------------------------------------------
    # 2. Pedicel Tip Mesh (x from -l_pedicel to 0)
    # -------------------------------------------------------------
    n_ped = int(n_x * 0.15)
    x_ped = np.linspace(-l_pedicel, 0, n_ped, endpoint=False)

    ped_grid_pts = []
    pedicel_tip_id = points.InsertNextPoint(-l_pedicel, 0.0, 0.0)

    for x in x_ped:
        eta = (x - (-l_pedicel)) / l_pedicel

        w_p = 0.2 + (w_top - 0.2) * (eta**1.2)
        t_p = 0.2 + (thickness - 0.2) * (eta**1.2)

        a_y = max(w_p / 2.0, 0.01)
        b_z = max(t_p / 2.0, 0.01)

        row = []
        for th in theta:
            cos_t = np.cos(th)
            sin_t = np.sin(th)
            y = a_y * np.sign(cos_t) * (np.abs(cos_t) ** (2.0 / 2.5))
            z = b_z * np.sign(sin_t) * (np.abs(sin_t) ** (2.0 / 2.5))
            pt_id = points.InsertNextPoint(x, y, z)
            row.append(pt_id)
        ped_grid_pts.append(row)

    all_grid_pts = ped_grid_pts + grid_pts
    total_x_layers = len(all_grid_pts)

    # -------------------------------------------------------------
    # 3. Assemble PolyData Cells & Normals
    # -------------------------------------------------------------
    for j in range(n_theta):
        j_next = (j + 1) % n_theta
        p1 = ped_grid_pts[0][j_next]
        p2 = ped_grid_pts[0][j]
        polys.InsertNextCell(3, [pedicel_tip_id, p1, p2])

    for i in range(total_x_layers - 1):
        for j in range(n_theta):
            j_next = (j + 1) % n_theta
            p1 = all_grid_pts[i][j]
            p2 = all_grid_pts[i][j_next]
            p3 = all_grid_pts[i + 1][j_next]
            p4 = all_grid_pts[i + 1][j]
            polys.InsertNextCell(4, [p1, p2, p3, p4])

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(polys)

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polydata)
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOn()
    normals.SplittingOff()
    normals.Update()

    return normals.GetOutput()


def visualize_maize_kernel_interactive(polydata):
    """Launches an interactive VTK 3D window."""
    # 1. Mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    mapper.ScalarVisibilityOff()

    # 2. Actor (Golden Corn Yellow material with Gouraud shading & specular gloss)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetColor(0.95, 0.75, 0.18)  # Golden yellow
    prop.SetSpecular(0.4)
    prop.SetSpecularPower(30.0)
    prop.SetAmbient(0.2)
    prop.SetDiffuse(0.8)
    prop.SetInterpolationToGouraud()

    # 3. Renderer & Gradient Background
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.12, 0.15, 0.22)
    renderer.SetBackground2(0.04, 0.06, 0.10)
    renderer.GradientBackgroundOn()

    # 4. Render Window
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    renWin.SetSize(1024, 768)
    renWin.SetWindowName("Figure 5b - Single Maize Kernel Model (Interactive VTK)")

    # 5. Interactor & Trackball Camera Control
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)

    # 6. Corner Coordinate Axes Orientation Marker
    axes_actor = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes_actor)
    axes_widget.SetInteractor(iren)
    axes_widget.SetViewport(0.0, 0.0, 0.2, 0.2)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # Position Camera
    renderer.ResetCamera()
    cam = renderer.GetActiveCamera()
    cam.Elevation(20)
    cam.Azimuth(30)

    print("Launching Interactive VTK Window...")
    print("Controls:")
    print("  - Left Mouse: Rotate")
    print("  - Middle Mouse / Shift+Left: Pan")
    print("  - Right Mouse / Scroll: Zoom")
    print("  - 'w' key: Wireframe mode")
    print("  - 's' key: Surface shading mode")
    print("  - 'r' key: Reset camera")
    print("  - 'e' / 'q' key: Exit window")

    renWin.Render()
    iren.Initialize()
    iren.Start()


if __name__ == "__main__":
    # Generate 3D PolyData mesh
    polydata = create_maize_kernel_mesh()
    print(f"Mesh generated: {polydata.GetNumberOfPoints()} vertices," f" {polydata.GetNumberOfCells()} faces.")

    # Export STL and VTP files
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetFileName("maize_kernel_fig5b.stl")
    stl_writer.SetInputData(polydata)
    stl_writer.Write()

    vtp_writer = vtk.vtkXMLPolyDataWriter()
    vtp_writer.SetFileName("maize_kernel_fig5b.vtp")
    vtp_writer.SetInputData(polydata)
    vtp_writer.Write()

    # Launch local interactive window
    visualize_maize_kernel_interactive(polydata)
