"""
Leaf animation example.

Creates JPEG frames first, assembles them into an MP4, removes the frame
folder to avoid clutter, and then opens the MP4 in VLC.
"""

import shutil
import subprocess
from pathlib import Path

import numpy as np
from leafs import example, long_leaf, mint, parsley

import plantbox as pb
from plantbox.visualisation.vtk_animate import AnimateRoots

DURATION_SECONDS = 8.0
TOTAL_SIM_DAYS = 30.0
FPS = 25
N_FRAMES = int(DURATION_SECONDS * FPS)


FRAME_DIR = Path("leafs_frames")
OUTPUT_MP4 = Path("leafs_animation.mp4")


def build_leaf_plant():

    plant = pb.Plant()

    seed_rp = pb.SeedRandomParameter(plant)
    seed_rp.subType = 0
    plant.setOrganRandomParameter(seed_rp)

    root_rp = pb.RootRandomParameter(plant)
    root_rp.subType = 1
    root_rp.lmax = 1
    root_rp.r = 0.1
    root_rp.theta = 0
    plant.setOrganRandomParameter(root_rp)

    stem_rp = pb.StemRandomParameter(plant)
    stem_rp.subType = 1
    stem_rp.la = 1.0
    stem_rp.lb = 5
    stem_rp.lmax = 7.5
    stem_rp.ln = 1
    stem_rp.theta = 0
    stem_rp.successor = [[2]]
    stem_rp.successorP = [[1]]
    stem_rp.successorOT = [[pb.leaf]]
    plant.setOrganRandomParameter(stem_rp)

    leaf_rp = pb.LeafRandomParameter(plant)
    leaf_rp.parametrisationType = 0
    leaf_rp.shapeType = 2
    leaf_rp.subType = 2
    leaf_rp.dx = 0.1
    leaf_rp.a = 0.05

    leaf_rp.tropismS = 0.0
    leaf_rp.theta = 90.0 / (2 * np.pi)

    # mint(leaf_rp)
    # parsley(leaf_rp)
    # example(leaf_rp)
    long_leaf(leaf_rp)
    plant.setOrganRandomParameter(leaf_rp)

    return plant


def play_in_vlc(video_path: Path):
    vlc = shutil.which("vlc")
    if vlc is None:
        print("VLC not found on PATH. Skipping playback.")
        return

    subprocess.Popen([vlc, str(video_path)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"Started VLC for {video_path}")


def ensure_ffmpeg():
    if shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg not found on PATH. Install ffmpeg to create MP4 files.")


def main():
    ensure_ffmpeg()

    # Start clean to avoid stale frames from previous runs.
    if FRAME_DIR.exists():
        shutil.rmtree(FRAME_DIR)

    plant = build_leaf_plant()
    plant.initialize(verbose=False)

    # Warm-up so camera bounds are not degenerate at frame 0.
    warmup_dt = 0.1
    elapsed = 0.0
    for _ in range(20):
        if plant.getNumberOfSegments() > 0:
            break
        plant.simulate(warmup_dt, verbose=False)
        elapsed += warmup_dt

    anim = AnimateRoots(plant)
    anim.root_name = "creationTime"
    anim.plant = True
    anim.avi_name = str(FRAME_DIR)
    anim.start(axis="v")

    dt = TOTAL_SIM_DAYS / N_FRAMES
    try:
        for _ in range(N_FRAMES):
            plant.simulate(dt, verbose=False)
            elapsed += dt
            anim.simtime = elapsed
            anim.update()

        anim.make_video(output_file=str(OUTPUT_MP4), fps=FPS)
    finally:
        if FRAME_DIR.exists():
            shutil.rmtree(FRAME_DIR)
            print(f"Removed temporary frame directory: {FRAME_DIR}")

    print(f"Created video: {OUTPUT_MP4.resolve()}")
    play_in_vlc(OUTPUT_MP4.resolve())


if __name__ == "__main__":
    main()
