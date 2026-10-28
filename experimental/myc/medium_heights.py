""" Estimates the medium heights in the split Petri dish of the experiment from the filled volumes,
    and where the opening in the acrylic frame lies relative to the medium surfaces.

    Setup of the experiment (split plate, 94 mm, Greiner Bio-One):
    - fungal compartment: 28 ml MSR, cellophane on top (the hyphal network grows on the cellophane)
    - root compartment: 25 ml 1%P MSR (5 ml + 20 ml), with the acrylic frame standing on the bottom against the central barrier
    - acrylic frame: 1 mm thick, trapezoid with top edge 88 mm, bottom edge 85.5 mm, height 12 mm,
      opening 50 x 2 mm located 2 mm below the top edge

    Not given in the text: inner diameter of the dish and thickness of the central barrier, so they are arguments.
    The simulation scripts call estimate_heights() with the diameter and barrier thickness of the model, then the
    modelled medium has the volume of the experiment. Run this file directly to see how the heights depend on them.
    All lengths in cm, volumes in ml = cm^3.
"""
import numpy as np

# experiment
V_fungal = 28.  # ml MSR in the fungal compartment
V_root = 25.  # ml 1%P MSR in the root compartment
frame_thickness = 0.1  # 1 mm
frame_bottom, frame_top, frame_height = 8.55, 8.8, 1.2  # trapezoid edges and height
opening_height = 0.2  # opening 50 x 2 mm ...
opening_below_top = 0.2  # ... located 2 mm below the top edge of the frame


def half_area(radius, wall):
    """ area of one compartment: disk of radius radius with x > wall/2 (central barrier of thickness wall at x = 0) """
    c = wall / 2
    return np.pi / 2 * radius * radius - (c * np.sqrt(radius * radius - c * c) + radius * radius * np.arcsin(c / radius))


def frame_volume_below(h, thickness, width_bottom, width_top, frame_height):
    """ volume of the acrylic frame below height h above the bottom (width changes linearly with height) """
    h = min(h, frame_height)
    return thickness * (width_bottom * h + (width_top - width_bottom) * h * h / (2 * frame_height))


def medium_height(volume, area, frame = None):
    """ height of a medium volume filled into a compartment of base area area, minus the immersed frame volume if given """
    if frame is None:
        return volume / area
    h = volume / area
    for i in range(50):  # fixed point iteration, the frame takes only a few percent of the volume
        h = (volume + frame_volume_below(h, *frame)) / area
    return h


def estimate_heights(diameter, wall):
    """ medium heights [cm] for a dish of (inner) diameter diameter with a central barrier of thickness wall:
        "fungal" and "root" medium height, "frame_volume" immersed in the root medium [ml], the opening
        "opening_bottom", "opening_top" above the bottom, and the compartment "area" [cm^2] """
    A = half_area(diameter / 2, wall)
    frame = (frame_thickness, frame_bottom, frame_top, frame_height)
    h_root = medium_height(V_root, A, frame)
    opening_top = frame_height - opening_below_top
    return {"fungal": medium_height(V_fungal, A), "root": h_root, "frame_volume": frame_volume_below(h_root, *frame),
            "opening_bottom": opening_top - opening_height, "opening_top": opening_top, "area": A}


def print_heights(h, diameter, wall):
    print(f"dish diameter {diameter * 10:.0f} mm, central barrier {wall * 10:.1f} mm: compartment area {h['area']:.1f} cm^2")
    print(f"   fungal medium {h['fungal'] * 10:.1f} mm, root medium {h['root'] * 10:.1f} mm (frame takes {h['frame_volume']:.2f} ml)")
    print(f"   opening {h['opening_bottom'] * 10:.0f} - {h['opening_top'] * 10:.0f} mm above the bottom: "
          f"{(h['opening_bottom'] - h['fungal']) * 10:+.1f} to {(h['opening_top'] - h['fungal']) * 10:+.1f} mm relative to the cellophane (fungal surface), "
          f"{(h['opening_bottom'] - h['root']) * 10:+.1f} to {(h['opening_top'] - h['root']) * 10:+.1f} mm relative to the root medium surface")


if __name__ == "__main__":
    wall = 0.1  # assumption, thickness of the central barrier of the split plate
    for diameter in [9.4, 9.0, 8.8]:  # 94 mm is the nominal size, the inner diameter is a bit smaller
        print_heights(estimate_heights(diameter, wall), diameter, wall)

    print("""
Conclusion for the model (makedishes):
- the medium is 7-9 mm high, not the 16 mm of the dish, the air above the medium is not part of the model: height = root medium height
  (the fungal surface is about 0.7 mm higher, which the model neglects by using one height)
- the opening lies at the medium surfaces, from about the cellophane level up to 1-2 mm above it:
  in the model it is at the top of the medium, z in [-opening_height, 0] (opening_at_top = True)
- the hyphae in the fungal compartment grow on the cellophane, i.e. in a thin layer at the top (hyphae_layer, e.g. 0.1 cm)
""")
