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

    # Allow importing test utilities from the repo (working dir or CI mounts)
    for root in (os.getcwd(), "/src", "/project"):
        tp = os.path.join(root, "test")
        if os.path.isdir(tp) and tp not in sys.path:
            sys.path.append(tp)

    from tools_image import compare_images_png, render_headless_png  # type: ignore

    # Diagnostics: data path and parameter file
    data_root = pb.data_path()
    param_dir = os.path.join(data_root, "structural", "plant")
    param_file = os.path.join(param_dir, "fspm2023.xml")
    print("data_root:", data_root)
    print("param_file exists:", os.path.exists(param_file))

    # Build deterministic plant
    plant = pb.Plant(1)  # type: ignore
    plant.readParameters(param_file)
    plant.initialize(False)
    plant.simulate(40, False)

    # Optional diagnostics
    try:
        ana = pb.SegmentAnalyser(plant)  # type: ignore
        print("nodes:", len(list(ana.nodes)), "segments:", len(list(ana.segments)))
    except Exception:
        pass

    # Render headless
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    out_png = "/tmp/generated_by_smoke.png"
    render_headless_png(
        plant, p_name="age", image_path=out_png, width=1000, height=int(1000 * 2.2), zoom=3.0
    )

    # Copy artifact into repo for inspection if available
    for dst_root in (os.getcwd(), "/src", "/project"):
        if not os.path.isdir(dst_root):
            continue
        dst = os.path.join(dst_root, "test", "golden", "generated_by_test.png")
        try:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with open(out_png, "rb") as fsrc, open(dst, "wb") as fdst:
                fdst.write(fsrc.read())
            print("wrote:", dst)
            break
        except Exception:
            pass

    # Compare with macOS golden reference
    root = next((p for p in (os.getcwd(), "/src", "/project") if os.path.isdir(p)), None)
    if not root:
        print("No repo root mount found")
        return 2
    mac_golden = os.path.join(root, "test", "golden", "macos", "example_plant_headless.png")
    ref_png = os.path.join(root, "test", "golden", "example_plant_headless.png")
    if not os.path.exists(mac_golden):
        print("Golden macOS image missing:", mac_golden)
        return 2
    try:
        os.makedirs(os.path.dirname(ref_png), exist_ok=True)
        with open(mac_golden, "rb") as fsrc, open(ref_png, "wb") as fdst:
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
