import os
import sys
import sysconfig


def main() -> int:
    # Prefer software rendering where possible
    os.environ.setdefault("LIBGL_ALWAYS_SOFTWARE", "1")
    os.environ.setdefault("VTK_USE_OFFSCREEN", "1")

    # Ensure we import the installed package, not the repo source
    site = sysconfig.get_paths()["platlib"]
    sys.path[:] = [site] + [
        p
        for p in sys.path
        if p != site and not p.startswith("/project") and not p.startswith("/src")
    ]

    import plantbox as pb  # installed from wheel

    # Allow importing test utilities from the repo (both possible mounts)
    test_paths = ["/src/test", "/project/test"]
    for tp in test_paths:
        if os.path.isdir(tp):
            sys.path.append(tp)
            break

    from tools_image import compare_images_png, render_headless_png  # type: ignore

    # Diagnostics: data path and parameter file
    data_root = pb.data_path()
    param_dir = os.path.join(data_root, "structural", "plant")
    param_file = os.path.join(param_dir, "fspm2023.xml")
    print("data_root:", data_root)
    print("param_file exists:", os.path.exists(param_file))

    # Build deterministic plant
    plant = pb.Plant(1)
    plant.readParameters(param_file)
    plant.initialize(False)
    plant.simulate(40, False)

    # Geometry diagnostics
    try:
        ana = pb.SegmentAnalyser(plant)
        num_nodes = len(list(ana.nodes))
        num_segments = len(list(ana.segments))
        print("nodes:", num_nodes, "segments:", num_segments)
    except Exception as e:  # pragma: no cover
        print("SegmentAnalyser error:", e)

    # Render headless with dimensions matching the golden
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    out_png = "/tmp/generated_by_smoke.png"
    width = 1000
    height = int(1000 * 2.2)
    render_headless_png(
        plant,
        p_name="age",
        image_path=out_png,
        width=width,
        height=height,
        zoom=3.0,
    )

    # Also copy image to mounted repo for inspection if available
    for dst in (
        "/src/test/golden/generated_by_test.png",
        "/project/test/golden/generated_by_test.png",
    ):
        try:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with open(out_png, "rb") as fsrc, open(dst, "wb") as fdst:
                fdst.write(fsrc.read())
            print("wrote:", dst)
            break
        except Exception:
            pass

    # Compare with golden: copy linux golden into generic path and compare (leave it in place for inspection)
    repo_roots = [p for p in ("/src", "/project") if os.path.isdir(p)]
    if not repo_roots:
        print("No repo root mount found; expected /src or /project")
        return 2
    root = repo_roots[0]
    linux_golden = os.path.join(root, "test", "golden", "linux", "example_plant_headless.png")
    ref_png = os.path.join(root, "test", "golden", "example_plant_headless.png")
    if not os.path.exists(linux_golden):
        print("Golden linux image missing:", linux_golden)
        return 2
    try:
        os.makedirs(os.path.dirname(ref_png), exist_ok=True)
        with open(linux_golden, "rb") as fsrc, open(ref_png, "wb") as fdst:
            fdst.write(fsrc.read())
        print("prepared reference:", ref_png)
    except Exception as e:
        print("Failed to prepare reference:", e)
        return 2

    sim = compare_images_png(out_png, ref_png)

    print("golden_similarity=", sim)
    return 0 if sim >= 0.9 else 2


if __name__ == "__main__":
    raise SystemExit(main())
