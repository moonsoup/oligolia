"""Generate Oligolia app icon — DNA helix, green on dark navy.

Produces:
  assets/icon.png    (1024x1024, for macOS / generic)
  assets/icon.ico    (multi-size, for Windows)
  assets/icon.icns   (macOS native, if iconutil available)
"""

import math
import os
from pathlib import Path

from PIL import Image, ImageDraw

SIZE = 1024
HALF = SIZE // 2
ASSETS = Path(__file__).parent


def draw_icon(size: int = 1024) -> Image.Image:
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    cx, cy = size // 2, size // 2
    r = int(size * 0.46)

    # Background circle — dark navy
    draw.ellipse([cx - r, cy - r, cx + r, cy + r], fill="#0a1628", outline="#1a3a5c", width=max(2, size // 120))

    # Draw double helix
    strand_r = size * 0.14       # radius of helix rotation
    strand_width = max(3, size // 100)
    rung_width = max(2, size // 140)
    helix_top = cy - int(size * 0.34)
    helix_bot = cy + int(size * 0.34)
    helix_height = helix_bot - helix_top
    turns = 2.5
    steps = 200
    GREEN = "#4ade80"
    BLUE = "#60a5fa"
    RUNG = "#94a3b8"

    def helix_point(t: float, phase_offset: float = 0) -> tuple[float, float]:
        angle = t * turns * 2 * math.pi + phase_offset
        x = cx + strand_r * math.sin(angle)
        y = helix_top + t * helix_height
        return x, y

    # Draw rungs first (behind strands)
    rung_count = int(turns * 6)
    for i in range(rung_count + 1):
        t = i / rung_count
        x1, y1 = helix_point(t, 0)
        x2, y2 = helix_point(t, math.pi)
        draw.line([(x1, y1), (x2, y2)], fill=RUNG, width=rung_width)

    # Strand 1 — green
    pts1 = [helix_point(i / steps, 0) for i in range(steps + 1)]
    for i in range(len(pts1) - 1):
        draw.line([pts1[i], pts1[i + 1]], fill=GREEN, width=strand_width)

    # Strand 2 — blue
    pts2 = [helix_point(i / steps, math.pi) for i in range(steps + 1)]
    for i in range(len(pts2) - 1):
        draw.line([pts2[i], pts2[i + 1]], fill=BLUE, width=strand_width)

    # Nucleotide dots on strand 1
    for i in range(0, steps + 1, steps // (rung_count)):
        x, y = helix_point(i / steps, 0)
        dot_r = max(4, size // 80)
        draw.ellipse([x - dot_r, y - dot_r, x + dot_r, y + dot_r], fill=GREEN)

    # Text "Og" watermark small in corner — skip for clean icon

    return img


def main() -> None:
    img = draw_icon(SIZE)
    png_path = ASSETS / "icon.png"
    img.save(png_path, "PNG")
    print(f"Saved {png_path}")

    # Windows ICO — multiple sizes
    ico_sizes = [(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
    ico_images = [img.resize(s, Image.LANCZOS) for s in ico_sizes]
    ico_path = ASSETS / "icon.ico"
    ico_images[0].save(ico_path, format="ICO", sizes=ico_sizes,
                       append_images=ico_images[1:])
    print(f"Saved {ico_path}")

    # macOS iconset (requires iconutil, available on macOS)
    iconset_dir = ASSETS / "icon.iconset"
    iconset_dir.mkdir(exist_ok=True)
    mac_sizes = [16, 32, 64, 128, 256, 512, 1024]
    for s in mac_sizes:
        resized = img.resize((s, s), Image.LANCZOS)
        resized.save(iconset_dir / f"icon_{s}x{s}.png", "PNG")
        if s <= 512:
            resized2 = img.resize((s * 2, s * 2), Image.LANCZOS)
            resized2.save(iconset_dir / f"icon_{s}x{s}@2x.png", "PNG")
    print(f"Saved iconset at {iconset_dir}")

    # Try to create .icns via iconutil (macOS only)
    import subprocess
    icns_path = ASSETS / "icon.icns"
    result = subprocess.run(
        ["iconutil", "-c", "icns", str(iconset_dir), "-o", str(icns_path)],
        capture_output=True,
    )
    if result.returncode == 0:
        print(f"Saved {icns_path}")
    else:
        print("iconutil not available — .icns skipped (OK on non-macOS)")


if __name__ == "__main__":
    main()
