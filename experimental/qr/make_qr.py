"""Generate QR codes that link to websites and export them as PNG files."""

import os

import matplotlib.pyplot as plt
import qrcode
from PIL import Image, ImageDraw, ImageFont

# Forschungszentrum Jülich corporate blue
FZJ_BLUE = (2, 61, 107)

_FONT_CANDIDATES = [
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
]


def _load_font(size: int) -> ImageFont.FreeTypeFont:
    """Load a bold TrueType font, falling back to Pillow's built-in font."""
    for path in _FONT_CANDIDATES:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default(size=size)


def _rounded_mask(size: tuple[int, int], radius: int) -> Image.Image:
    """Create a rounded-rectangle alpha mask."""
    mask = Image.new("L", size, 0)
    ImageDraw.Draw(mask).rounded_rectangle([(0, 0), (size[0] - 1, size[1] - 1)], radius=radius, fill=255)
    return mask


def make_qr(name: str, hyperlink: str, outdir: str = "output", show: bool = False) -> str:
    """Create a styled QR code for a hyperlink and save it as a PNG file.

    The QR code uses the FZ Jülich corporate blue, sits inside a rounded
    box of the same colour, and shows ``name`` as a white caption below.

    Args:
        name: Name of the output file (without extension), shown as caption.
        hyperlink: The URL the QR code should point to.
        outdir: Directory where the PNG is saved (created if missing).
        show: If True, display the result in a matplotlib window.

    Returns:
        The path of the saved PNG file.
    """
    qr = qrcode.QRCode(
        version=None,  # auto-size
        error_correction=qrcode.constants.ERROR_CORRECT_M,
        box_size=10,
        border=4,
    )
    qr.add_data(hyperlink)
    qr.make(fit=True)
    qr_img = qr.make_image(fill_color=FZJ_BLUE, back_color="white").convert("RGBA")

    padding = qr_img.width // 10
    padding_x = qr_img.width // 10  # narrower margins left and right of the box
    margin_top = qr_img.width // 40  # extra blue margin at the top of the box
    radius = qr_img.width // 8
    caption_space = qr_img.width // 4

    box_w = qr_img.width + 2 * padding_x
    box_h = qr_img.height + 2 * padding + caption_space + margin_top

    # rounded box with FZJ blue background (transparent corners outside)
    img = Image.new("RGBA", (box_w, box_h), (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle([(0, 0), (box_w - 1, box_h - 1)], radius=radius, fill=FZJ_BLUE + (255,))

    # paste the QR code with slightly rounded corners
    img.paste(qr_img, (padding_x, padding + margin_top), _rounded_mask(qr_img.size, radius // 3))

    # white caption below the QR code, scaled to fit the box width
    font_size = qr_img.width // 9
    font = _load_font(font_size)
    bbox = draw.textbbox((0, 0), name, font=font)
    text_w = bbox[2] - bbox[0]
    while text_w > box_w - padding and font_size > 8:
        font_size -= 2
        font = _load_font(font_size)
        bbox = draw.textbbox((0, 0), name, font=font)
        text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    text_x = (box_w - text_w) / 2 - bbox[0]
    text_y = margin_top + padding + qr_img.height + (caption_space - text_h) / 2 - bbox[1]
    draw.text((text_x, text_y), name, font=font, fill="white")

    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{name}.png")
    img.save(path)

    if show:
        plt.imshow(img)
        plt.axis("off")
        plt.show()

    return path


if __name__ == "__main__":

    make_qr("CPlantBox", "https://www.cplantbox.com", show=False)
    make_qr("CPlantBox WebApp", "https://cplantbox.fz-juelich.de/", show=False)
    make_qr("CPlantBox GitHub", "https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox", show=False)
    make_qr("dumux-rosi", "https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi", show=False)
