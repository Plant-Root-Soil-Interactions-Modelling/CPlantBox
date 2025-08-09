import os

import plantbox as _pb  # type: ignore
from test.tools_image import compare_images_png, render_headless_png


def test_golden_headless_render(tmp_path):
    # Paths
    project_root = os.path.dirname(os.path.dirname(__file__))
    golden_dir = os.path.join(project_root, "test", "golden")
    golden_png = os.path.join(golden_dir, "example_plant_headless.png")

    # Prepare output path (temporary)
    # out_png = tmp_path / "render.png"
    out_png = os.path.join(golden_dir, "generated_by_test.png")

    # Build deterministic plant
    pb = _pb  # type: ignore
    plant = getattr(pb, "Plant")(1)
    path = os.path.join(pb.data_path(), "structural", "plant") + "/"
    name = "fspm2023"
    plant.readParameters(path + name + ".xml")
    plant.initialize(False)
    plant.simulate(40, False)

    # Render headless
    render_headless_png(plant, p_name="age", image_path=str(out_png))

    # Compare with golden
    assert os.path.exists(golden_png), "Golden image not found; generate via tut_1_3_headless.py"
    sim = compare_images_png(str(out_png), golden_png)

    # Be generous to account for platform differences; adjust as needed
    assert sim >= 0.9987274, f"Headless render deviates from golden (similarity={sim:.4f})"
