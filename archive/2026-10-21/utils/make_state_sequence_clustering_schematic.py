from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path("/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi")
OUT = ROOT / "output" / "state_sequence_clustering_schematic_20260818.png"

FONT = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"


def get_font(size: int, bold: bool = False):
    return ImageFont.truetype(FONT_BOLD if bold else FONT, size=size)


def rounded(draw, xy, fill, outline, radius=24, width=2):
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=width)


def wrapped(draw, text, xy, max_width, font, fill, line_gap=5):
    words = text.split()
    lines = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if draw.textlength(trial, font=font) <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)

    x, y = xy
    box = draw.textbbox((0, 0), "Ag", font=font)
    line_h = box[3] - box[1]
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        y += line_h + line_gap
    return y


def arrow(draw, start, end, color, width=5, head=12):
    draw.line((start, end), fill=color, width=width)
    x1, y1 = start
    x2, y2 = end
    if x2 == x1:
        return
    draw.polygon(
        [
            (x2, y2),
            (x2 - head, y2 - head // 2),
            (x2 - head, y2 + head // 2),
        ],
        fill=color,
    )


def elbow_arrow(draw, start, end, mid_y, color, width=5, head=12):
    x1, y1 = start
    x2, y2 = end
    draw.line((x1, y1, x1, mid_y), fill=color, width=width)
    draw.line((x1, mid_y, x2, mid_y), fill=color, width=width)
    draw.line((x2, mid_y, x2, y2), fill=color, width=width)
    draw.polygon(
        [
            (x2, y2),
            (x2 - head // 2, y2 - head),
            (x2 + head // 2, y2 - head),
        ],
        fill=color,
    )


def draw_sequence(draw, x, y, block_w, block_h, colors, labels=None, outline="#ffffff"):
    for i, color in enumerate(colors):
        x0 = x + i * block_w
        rounded(draw, (x0, y, x0 + block_w - 6, y + block_h), color, outline, radius=8, width=1)
        if labels:
            bbox = draw.textbbox((0, 0), labels[i], font=get_font(16, True))
            tw = bbox[2] - bbox[0]
            th = bbox[3] - bbox[1]
            draw.text((x0 + (block_w - 6 - tw) / 2, y + (block_h - th) / 2 - 1), labels[i], font=get_font(16, True), fill="#17324d")


def main():
    W, H = 2400, 1400
    img = Image.new("RGB", (W, H), "#f6f8fb")
    draw = ImageDraw.Draw(img)

    navy = "#17324d"
    blue = "#2f6fb2"
    slate = "#4d6077"
    pale = "#d8e1ee"
    white = "#ffffff"
    green = "#9dbb6c"      # room air
    yellow = "#f5d36f"     # low-flow oxygen
    teal = "#86b9c6"       # NIV
    red = "#d86a74"        # IMV
    light_panel = "#fbfcfe"
    accent = "#6a7e95"

    section_font = get_font(28, True)
    body_font = get_font(22, False)
    small_font = get_font(18, False)
    tiny_font = get_font(16, False)

    # legend in its own box
    legend_box = (1250, 118, 2225, 218)
    rounded(draw, legend_box, light_panel, pale, radius=18, width=2)
    draw.text((1275, 130), "Respiratory support states", font=small_font, fill=navy)
    legend_items = [
        (green, "Room air", 1275, 166),
        (yellow, "Low-flow oxygen", 1455, 166),
        (teal, "Noninvasive ventilation", 1690, 166),
        (red, "Invasive ventilation", 1958, 166),
    ]
    for color, label, lx, ly in legend_items:
        rounded(draw, (lx, ly - 2, lx + 24, ly + 22), color, color, radius=6, width=1)
        draw.text((lx + 34, ly), label, font=tiny_font, fill=slate)

    # panels
    panels = {
        "step1": (50, 250, 760, 700),
        "step2": (825, 250, 1535, 700),
        "step3": (1600, 250, 2350, 700),
        "step4": (360, 805, 1190, 1340),
        "step5": (1310, 805, 2140, 1340),
    }
    for box in panels.values():
        rounded(draw, box, light_panel, pale, radius=28, width=3)

    # arrows between panels
    arrow(draw, (760, 475), (825, 475), accent)
    arrow(draw, (1535, 475), (1600, 475), accent)
    elbow_arrow(draw, (1975, 700), (775, 805), 750, accent)
    arrow(draw, (1190, 1073), (1310, 1073), accent)

    # Step 1
    x1, y1, x2, y2 = panels["step1"]
    draw.text((x1 + 20, y1 + 18), "1. Encode each patient as an", font=section_font, fill=navy)
    draw.text((x1 + 20, y1 + 54), "ordered state sequence", font=section_font, fill=navy)
    wrapped(
        draw,
        "For each hospitalization, the highest respiratory support level in each hour becomes one state in a fixed-length sequence.",
        (x1 + 25, y1 + 98),
        620,
        body_font,
        slate,
    )
    draw.text((x1 + 25, y1 + 172), "Hours from ICU admission", font=small_font, fill=slate)
    hour_labels = [str(i) for i in range(1, 9)]
    draw_sequence(draw, x1 + 245, y1 + 154, 46, 34, [white] * 8, labels=hour_labels, outline=pale)
    seqs = [
        ("Patient A", [green, green, yellow, yellow, teal, teal, yellow, green]),
        ("Patient B", [yellow, yellow, teal, teal, teal, yellow, yellow, green]),
        ("Patient C", [red, red, red, red, red, yellow, yellow, teal]),
    ]
    sy = y1 + 222
    for name, seq in seqs:
        draw.text((x1 + 30, sy + 8), name, font=body_font, fill=navy)
        draw_sequence(draw, x1 + 245, sy, 46, 42, seq)
        sy += 74

    # Step 2
    x1, y1, x2, y2 = panels["step2"]
    draw.text((x1 + 25, y1 + 20), "2. Compare every pair of sequences", font=section_font, fill=navy)
    wrapped(
        draw,
        "Optimal matching measures how much work it takes to transform one sequence into another using substitutions, insertions, and deletions.",
        (x1 + 25, y1 + 62),
        640,
        body_font,
        slate,
    )
    # two sequences and alignment
    draw.text((x1 + 45, y1 + 155), "Example comparison", font=small_font, fill=slate)
    seq_a = [green, green, yellow, teal, teal, yellow, green]
    seq_b = [yellow, yellow, teal, teal, yellow, yellow, green]
    draw.text((x1 + 45, y1 + 205), "A", font=body_font, fill=navy)
    draw_sequence(draw, x1 + 90, y1 + 190, 55, 40, seq_a)
    draw.text((x1 + 45, y1 + 275), "B", font=body_font, fill=navy)
    draw_sequence(draw, x1 + 90, y1 + 260, 55, 40, seq_b)
    # highlight differences
    diff_x = [2, 0, 5]
    for idx in diff_x:
        px = x1 + 90 + idx * 55
        draw.rectangle((px - 3, y1 + 187, px + 49, y1 + 303), outline=blue, width=3)
    draw.text((x1 + 90, y1 + 338), "Similar transitions -> lower cost", font=small_font, fill=slate)
    draw.text((x1 + 90, y1 + 370), "More mismatches or timing shifts -> higher cost", font=small_font, fill=slate)

    # Step 3
    x1, y1, x2, y2 = panels["step3"]
    draw.text((x1 + 25, y1 + 20), "3. Build a pairwise distance matrix", font=section_font, fill=navy)
    wrapped(
        draw,
        "Distances summarize how dissimilar every patient is from every other patient after sequence alignment.",
        (x1 + 25, y1 + 62),
        640,
        body_font,
        slate,
    )
    grid_x = x1 + 135
    grid_y = y1 + 178
    cell = 54
    labels = ["A", "B", "C", "D", "E"]
    vals = [
        [0, 2, 7, 3, 8],
        [2, 0, 6, 2, 7],
        [7, 6, 0, 5, 2],
        [3, 2, 5, 0, 6],
        [8, 7, 2, 6, 0],
    ]
    for i, lab in enumerate(labels):
        draw.text((grid_x + 56 + i * cell + 15, grid_y - 36), lab, font=body_font, fill=navy)
        draw.text((grid_x + 12, grid_y + 12 + i * cell), lab, font=body_font, fill=navy)
    for i in range(5):
        for j in range(5):
            val = vals[i][j]
            if val <= 2:
                fill = "#dff1e3"
            elif val <= 4:
                fill = "#eef4fb"
            elif val <= 6:
                fill = "#fdeecf"
            else:
                fill = "#f7d9dd"
            rounded(draw, (grid_x + 48 + j * cell, grid_y + i * cell, grid_x + 48 + (j + 1) * cell - 6, grid_y + (i + 1) * cell - 6), fill, white, radius=8, width=1)
            draw.text((grid_x + 48 + j * cell + 17, grid_y + i * cell + 11), str(val), font=small_font, fill=navy)
    # Step 4
    x1, y1, x2, y2 = panels["step4"]
    draw.text((x1 + 25, y1 + 20), "4. Cluster patients with similar sequences", font=section_font, fill=navy)
    wrapped(
        draw,
        "Hierarchical agglomerative clustering groups patients whose aligned trajectories have small pairwise distances.",
        (x1 + 25, y1 + 62),
        700,
        body_font,
        slate,
    )
    # simple dendrogram
    base_x = x1 + 85
    base_y = y1 + 210
    row_h = 58
    leaves = ["A", "B", "D", "C", "E"]
    for i, lab in enumerate(leaves):
        y = base_y + i * row_h
        draw.text((base_x, y - 10), lab, font=body_font, fill=navy)
        draw.line((base_x + 35, y, base_x + 140, y), fill=accent, width=4)
    # merge A-B
    draw.line((base_x + 140, base_y, base_x + 140, base_y + row_h), fill=accent, width=4)
    draw.line((base_x + 140, base_y + row_h / 2, base_x + 245, base_y + row_h / 2), fill=accent, width=4)
    # merge AB-D
    y_d = base_y + 2 * row_h
    draw.line((base_x + 140, y_d, base_x + 245, y_d), fill=accent, width=4)
    draw.line((base_x + 245, base_y + row_h / 2, base_x + 245, y_d), fill=accent, width=4)
    draw.line((base_x + 245, base_y + 1.25 * row_h, base_x + 375, base_y + 1.25 * row_h), fill=accent, width=4)
    # merge C-E
    y_c = base_y + 3 * row_h
    y_e = base_y + 4 * row_h
    draw.line((base_x + 140, y_c, base_x + 140, y_e), fill=accent, width=4)
    draw.line((base_x + 140, (y_c + y_e) / 2, base_x + 295, (y_c + y_e) / 2), fill=accent, width=4)
    # final merge
    draw.line((base_x + 375, base_y + 1.25 * row_h, base_x + 460, base_y + 1.25 * row_h), fill=accent, width=4)
    draw.line((base_x + 295, (y_c + y_e) / 2, base_x + 460, (y_c + y_e) / 2), fill=accent, width=4)
    draw.line((base_x + 460, base_y + 1.25 * row_h, base_x + 460, (y_c + y_e) / 2), fill=accent, width=4)
    draw.line((x1 + 585, y1 + 245, x1 + 585, y1 + 470), fill=blue, width=4)
    draw.text((x1 + 600, y1 + 345), "k = 2", font=body_font, fill=blue)

    # Step 5
    x1, y1, x2, y2 = panels["step5"]
    draw.text((x1 + 25, y1 + 20), "5. Interpret the resulting phenotypes", font=section_font, fill=navy)
    wrapped(
        draw,
        "Each cluster represents a recurring longitudinal pattern of respiratory support rather than a single time-point severity score.",
        (x1 + 25, y1 + 62),
        700,
        body_font,
        slate,
    )
    clusters = [
        ("Cluster 1", "Resolved trajectory", [green, green, yellow, yellow, green, green, green], "#e7f3ea"),
        ("Cluster 2", "Persistent high-acuity trajectory", [red, red, red, red, yellow, yellow, teal], "#faebec"),
    ]
    cy = y1 + 165
    for title, subtitle, seq, fill in clusters:
        rounded(draw, (x1 + 45, cy, x2 - 45, cy + 125), fill, pale, radius=20, width=2)
        draw.text((x1 + 70, cy + 18), title, font=body_font, fill=navy)
        draw.text((x1 + 70, cy + 52), subtitle, font=small_font, fill=slate)
        draw_sequence(draw, x1 + 370, cy + 34, 48, 38, seq)
        cy += 155

    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT, quality=95)
    print(OUT)


if __name__ == "__main__":
    main()
