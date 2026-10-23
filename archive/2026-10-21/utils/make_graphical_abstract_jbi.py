from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path("/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi")
OUT = ROOT / "output" / "jbi_graphical_abstract_20260812.png"
CONSORT = ROOT / "output" / "pooled_20260430" / "figure_pooled_consort_diagram.png"
TRAJ = ROOT / "output" / "pooled_20260430" / "figure_pooled_resp_support_trajectories.png"
OUTCOMES = ROOT / "output" / "pooled_20260430" / "figure_phenotype_binary_outcomes.png"

FONT = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"


def get_font(size: int, bold: bool = False):
    return ImageFont.truetype(FONT_BOLD if bold else FONT, size=size)


def rounded(draw, xy, fill, outline, radius=24, width=2):
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=width)


def add_wrapped(draw, text, xy, max_width, font, fill, line_gap=6):
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
    line_box = draw.textbbox((0, 0), "Ag", font=font)
    line_h = line_box[3] - line_box[1]
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        y += line_h + line_gap
    return y


def place_image(canvas, path, target_box, pad=0, background=None, crop=None):
    img = Image.open(path).convert("RGB")
    if crop:
        img = img.crop(crop)
    x1, y1, x2, y2 = target_box
    if background:
        bg = Image.new("RGB", (x2 - x1, y2 - y1), background)
        canvas.paste(bg, (x1, y1))
    tw = (x2 - x1) - 2 * pad
    th = (y2 - y1) - 2 * pad
    ratio = min(tw / img.width, th / img.height)
    nw = int(img.width * ratio)
    nh = int(img.height * ratio)
    img = img.resize((nw, nh))
    px = x1 + pad + (tw - nw) // 2
    py = y1 + pad + (th - nh) // 2
    canvas.paste(img, (px, py))


def main():
    W, H = 2800, 1700
    img = Image.new("RGB", (W, H), "#f5f7fb")
    draw = ImageDraw.Draw(img)

    navy = "#17324d"
    blue = "#2f6fb2"
    slate = "#4f6078"
    pale = "#d6e0ee"
    white = "#ffffff"
    green = "#2f9657"
    mild_blue = "#3875bb"
    gold = "#d7a400"
    orange = "#cb622a"

    title_font = get_font(56, True)
    sub_font = get_font(28, False)
    section_font = get_font(30, True)
    body_font = get_font(25, False)
    small_font = get_font(22, False)
    label_font = get_font(24, True)
    stat_font = get_font(30, True)

    draw.text((70, 40), "Graphical Abstract", font=section_font, fill=blue)
    draw.text(
        (70, 88),
        "Federated temporal phenotyping of respiratory trajectories in intensive care patients with lung cancer",
        font=title_font,
        fill=navy,
    )
    draw.text(
        (70, 156),
        "Journal of Biomedical Informatics framing: multicenter harmonized ICU EHR data, temporal sequence modeling, and clinically interpretable federated phenotypes",
        font=sub_font,
        fill=slate,
    )

    intro = (60, 230, 760, 1045)
    methods = (800, 230, 1360, 1045)
    results = (1400, 230, 2740, 1495)
    discussion = (60, 1095, 1360, 1495)

    for panel in (intro, methods, results, discussion):
        rounded(draw, panel, white, pale, radius=30, width=3)

    # Intro
    draw.text((90, 265), "Introduction and Study Cohort", font=section_font, fill=navy)
    add_wrapped(
        draw,
        "Patients with lung cancer admitted to the ICU have highly heterogeneous courses of acute respiratory failure. Most prognostic approaches reduce ICU data to static baseline severity scores and do not capture the time-evolving respiratory course.",
        (90, 320),
        620,
        body_font,
        slate,
        line_gap=6,
    )
    rounded(draw, (88, 500, 732, 995), "#f8fafc", pale, radius=24, width=2)
    place_image(
        img,
        CONSORT,
        (105, 525, 715, 920),
        pad=8,
        background="#f8fafc",
        crop=(0, 430, 4160, 2580),
    )
    rounded(draw, (100, 932, 260, 985), "#eef4fb", blue, radius=18, width=2)
    draw.text((126, 946), "Study cohort", font=small_font, fill=navy)
    draw.text((285, 941), "11,102 ICU hospitalizations", font=label_font, fill=navy)
    draw.text((285, 973), "10 U.S. health systems in CLIF", font=small_font, fill=slate)

    # Methods
    draw.text((830, 265), "Methods", font=section_font, fill=navy)
    steps = [
        ("1", "Federated CLIF data model", "Longitudinal ICU EHR data were harmonized locally across 10 health systems."),
        ("2", "72-hour temporal representation", "Hourly respiratory support and acute respiratory failure states were encoded from ICU admission."),
        ("3", "State-sequence clustering", "Optimal matching plus hierarchical clustering identified consensus respiratory phenotypes."),
        ("4", "Validation", "Phenotypes were tested against mortality, ICU length of stay, and post-72-hour outcomes."),
    ]
    y = 345
    for num, head, desc in steps:
        rounded(draw, (835, y, 888, y + 53), blue, blue, radius=18, width=1)
        draw.text((855, y + 9), num, font=small_font, fill="white")
        draw.text((910, y + 2), head, font=label_font, fill=navy)
        y = add_wrapped(draw, desc, (910, y + 34), 390, body_font, slate, line_gap=5) + 26
        if num != "4":
            draw.line((861, y - 10, 861, y + 14), fill=pale, width=5)

    rounded(draw, (835, 850, 1325, 995), "#eef4fb", pale, radius=22, width=2)
    draw.text((860, 882), "Informatics novelty", font=label_font, fill=navy)
    add_wrapped(
        draw,
        "A privacy-preserving temporal phenotyping framework generated reproducible, clinically interpretable respiratory subtypes from harmonized multicenter ICU data.",
        (860, 922),
        430,
        body_font,
        slate,
        line_gap=5,
    )

    # Results
    draw.text((1430, 265), "Results", font=section_font, fill=navy)
    draw.text((1440, 315), "Four respiratory trajectory phenotypes", font=label_font, fill=navy)
    rounded(draw, (1435, 350, 2705, 1020), "#fbfcfe", pale, radius=24, width=2)
    place_image(img, TRAJ, (1450, 365, 2690, 1005), pad=8, background="#fbfcfe")

    pill_y = 1045
    pills = [
        ("Resolved", "n=4,247", green, "#eaf5ee"),
        ("Mild", "n=3,908", mild_blue, "#eaf1fb"),
        ("Moderate", "n=1,332", gold, "#fff7df"),
        ("Persistent", "n=1,615", orange, "#fbede6"),
    ]
    x = 1445
    for label, n, color, fill in pills:
        rounded(draw, (x, pill_y, x + 270, pill_y + 95), fill, color, radius=24, width=3)
        draw.text((x + 18, pill_y + 18), label, font=label_font, fill=color)
        draw.text((x + 18, pill_y + 54), n, font=small_font, fill=slate)
        x += 305

    draw.text((1440, 1175), "Outcome gradient", font=label_font, fill=navy)
    rounded(draw, (1435, 1210, 2705, 1475), "#fbfcfe", pale, radius=24, width=2)
    place_image(
        img,
        OUTCOMES,
        (1450, 1225, 2320, 1455),
        pad=4,
        background="#fbfcfe",
        crop=(260, 0, 3800, 1980),
    )

    stats = [
        ("Resolved", "13%", green),
        ("Mild", "30%", mild_blue),
        ("Moderate", "60%", gold),
        ("Persistent", "66%", orange),
    ]
    sy = 1240
    for label, value, color in stats:
        rounded(draw, (2350, sy, 2670, sy + 40), white, pale, radius=14, width=1)
        draw.ellipse((2362, sy + 7, 2388, sy + 33), fill=color)
        draw.text((2402, sy + 4), label, font=small_font, fill=navy)
        draw.text((2588, sy + 3), value, font=small_font, fill=color)
        sy += 48

    # Discussion
    draw.text((90, 1105), "Discussion", font=section_font, fill=navy)
    rounded(draw, (88, 1160, 1330, 1360), "#eef4fb", blue, radius=24, width=2)
    add_wrapped(
        draw,
        "Temporal phenotyping revealed four clinically meaningful respiratory subtypes and identified a Moderate phenotype with high subsequent mortality that was not readily isolated by baseline SOFA or Charlson scores alone.",
        (120, 1200),
        1160,
        body_font,
        navy,
        line_gap=6,
    )
    draw.line((70, 1560, 2730, 1560), fill=pale, width=2)
    add_wrapped(
        draw,
        "Take-home message: A harmonized longitudinal ICU data model enabled federated sequence-based phenotyping of acute respiratory failure in lung cancer and produced interpretable subtypes with strong prognostic signal.",
        (90, 1585),
        2560,
        small_font,
        slate,
        line_gap=5,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    img.save(OUT, quality=95)
    print(OUT)


if __name__ == "__main__":
    main()
