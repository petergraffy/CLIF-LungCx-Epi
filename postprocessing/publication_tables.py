import argparse
import csv
from pathlib import Path
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

parser = argparse.ArgumentParser(description="Render manuscript Tables 1 and 2 from aggregate estimates.")
parser.add_argument("pool", type=Path)
args = parser.parse_args()
source = args.pool
dest = source / "tables"
dest.mkdir(exist_ok=True)

cohort = 'primary_phenotypes'
groups = ['Overall', 'P1', 'P2', 'P3', 'P4']
names = dict(zip(groups, ['Overall', 'P1 Shift toward\nroom air', 'P2 Sustained\nlow-flow oxygen', 'P3 High-flow\npredominant', 'P4 Sustained\ninvasive ventilation']))

def load(path):
    with (source / path).open() as f:
        return list(csv.DictReader(f))

denoms = load('publication_tables/denominators_by_site.csv')
ns = {g: sum(int(r['n']) for r in denoms if r['cohort'] == cohort and r['group'] == g) for g in groups + ['I']}
assert ns['Overall'] == sum(ns[g] for g in groups[1:] + ['I'])
cats = {(r['group'],r['variable'],r['category']):r for r in load('publication_tables/pooled_categorical.csv') if r['cohort'] == cohort}
cont = {(r['group'],r['variable']):r for r in load('publication_tables/pooled_continuous.csv') if r['cohort'] == cohort}

def build(title, headers, widths, rows, notes, landscape=True):
    doc = Document()
    sec = doc.sections[0]
    if landscape:
        sec.orientation = WD_ORIENT.LANDSCAPE
        sec.page_width, sec.page_height = Inches(11), Inches(8.5)
    sec.top_margin, sec.bottom_margin = Inches(.52), Inches(.52)
    sec.left_margin, sec.right_margin = Inches(.55), Inches(.55)
    style = doc.styles['Normal']
    style.font.name, style.font.size = 'Arial', Pt(9)
    style.paragraph_format.space_after = Pt(3)
    title_style = doc.styles['Title']
    title_style.font.name, title_style.font.size = 'Arial', Pt(15)
    title_style.font.bold = True
    title_style.font.color.rgb = RGBColor(0,0,0)
    title_style.paragraph_format.space_after = Pt(8)
    prop = title_style.element.get_or_add_pPr()
    for border in prop.findall(qn("w:pBdr")):
        prop.remove(border)
    doc.add_paragraph(title, 'Title')
    t = doc.add_table(rows=0, cols=len(headers))
    t.autofit = False
    for col, w in zip(t.columns,widths):
        col.width = Inches(w)
    for values, kind in [(headers,'header'), *rows]:
        row = t.add_row()
        prop = row._tr.get_or_add_trPr()
        prop.append(OxmlElement('w:cantSplit'))
        if kind == 'header':
            prop.append(OxmlElement('w:tblHeader'))
        if kind == 'section':
            row.cells[0].merge(row.cells[-1])
            cells = [row.cells[0]]
            used_widths = [sum(widths)]
        else:
            cells, used_widths = row.cells, widths
        for i,(cell,value,width) in enumerate(zip(cells, values, used_widths)):
            cell.width = Inches(width)
            if kind in ('header','section'):
                shd = OxmlElement('w:shd')
                shd.set(qn('w:fill'), 'E6EDF3' if kind == 'header' else 'F1F4F7')
                cell._tc.get_or_add_tcPr().append(shd)
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT if i == 0 else WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.line_spacing = 1.0
            p.add_run(value).bold = kind in ('header','section')
    for i,note in enumerate(notes):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(7 if i==0 else 0)
        p.add_run(note).font.size = Pt(8.5)
    return doc

def categorical(label, var, *categories):
    values=[]
    for g in groups:
        selected=[cats[(g,var,cat)] for cat in categories if (g,var,cat) in cats]
        assert selected and all(int(r['denominator_n']) == ns[g] for r in selected), (var,g,categories)
        n=sum(int(r['n']) for r in selected)
        values.append(f'{n:,} ({100*n/ns[g]:.1f})')
    return ([label,*values],'body')

def continuous(label, var, observed=False):
    values=[]
    for g in groups:
        r=cont[(g,var)]
        assert int(r['denominator_n']) == ns[g] and r['mean'] and r['sd']
        stat = f"{float(r['mean']):.1f} ({float(r['sd']):.1f})"
        if observed:
            stat += f"\nn = {int(r['observed_n']):,}"
        values.append(stat)
    return ([label,*values],'body')

def section(label):
    return ([label],'section')

rows = [
    section('Demographics and clinical characteristics'),
    continuous('Age, years', 'age_years'),
    categorical('Male sex', 'sex_group','Male'),
    categorical('White race', 'race_group','White'),
    categorical('Black race', 'race_group','Black'),
    categorical('Other or unknown race', 'race_group','Other known race','Unknown'),
    continuous('Charlson score excluding malignancy','charlson_noncancer'),
    categorical('Advanced cancer proxy','advanced_cancer','1'),
    categorical('Full code at ICU entry','code_at_admission','full'),
    categorical('Documented limitation at ICU entry','code_at_admission','dnr_without_documented_dni','dni','comfort'),
    categorical('Unknown code status at ICU entry','code_at_admission','unknown'),
    continuous('Modified SOFA','sofa_total',True),
    continuous('Adapted APACHE II','apache_adapted_total',True),
    section('First 72 hours after ICU entry'),
    categorical('Any high-flow nasal cannula','any_hfnc','TRUE'),
    categorical('Any noninvasive positive-pressure ventilation','any_nippv','TRUE'),
    categorical('Any invasive mechanical ventilation','any_imv','1'),
    categorical('Any documented ARF evidence','measured_arf','TRUE'),
    section('In-hospital outcomes through day 28 after landmark'),
    categorical('Death','outcome28','death'),
    categorical('Hospice discharge','outcome28','hospice'),
    categorical('Death or hospice discharge','outcome28','death','hospice'),
    categorical('Live discharge without hospice','outcome28','live_discharge'),
]
t1 = build('Table 1. Characteristics of landmark hospitalizations by respiratory phenotype',
    ['Characteristic',*(f'{names[g]}\nN = {ns[g]:,}' for g in groups)], [3.15,1.35,1.35,1.35,1.35,1.35], rows,
    [f"Values are n (%) or mean (SD). Percentages use the column total. The n shown under each severity score is the number with an available score. The overall column includes {ns['I']:,} hospitalizations with insufficient respiratory support information.",
     'ARF = acute respiratory failure. HFNC = high-flow nasal cannula. SOFA = Sequential Organ Failure Assessment. APACHE = Acute Physiology and Chronic Health Evaluation.'])
t1.save(dest / 'Table 1 - Cohort Characteristics.docx')

endpoints = {'death':'In-hospital death','hospice':'Hospice discharge','death_or_hospice':'Death or hospice discharge','live_discharge':'Live discharge without hospice'}
assocs = {(r['endpoint'],r['phenotype']): r for r in load('pooled_primary.csv') if r['model']=='core'}
rows=[]
for endpoint, label in endpoints.items():
    for idx,g in enumerate(groups[2:]):
        r=assocs[(endpoint,g)]
        assert r['reference']=='P1' and r['status']=='ok'
        hr=f"{float(r['hazard_ratio']):.2f} ({float(r['hr_low']):.2f}–{float(r['hr_high']):.2f})"
        pi=f"{float(r['hr_prediction_low']):.2f}–{float(r['hr_prediction_high']):.2f}"
        rows.append(([label if idx==0 else '', g, r['k'], hr, f"{float(r['i2_percent']):.1f}", pi],'body'))
t2=build('Table 2. Adjusted associations between respiratory phenotype and hospital outcomes',
    ['Outcome','Phenotype','Sites','Adjusted HR (95% CI)','I², %','95% prediction interval'], [2.55,.95,.75,2.15,.85,2.65], rows,
    ['P1 (shift toward room air) is the reference for every comparison. P2 = sustained low-flow oxygen. P3 = high-flow predominant. P4 = sustained invasive ventilation.',
     'Site-specific cause-specific Cox models used follow-up through 28 days after the 72-hour landmark. Models adjusted for age, sex, race, calendar year, Charlson score excluding malignancy, advanced cancer, and code status at ICU entry.',
     'Estimates were pooled using restricted maximum likelihood random-effects meta-analysis with modified Hartung–Knapp confidence intervals. HR = hazard ratio. CI = confidence interval.'])
t2.save(dest / 'Table 2 - Adjusted Cox Models.docx')
