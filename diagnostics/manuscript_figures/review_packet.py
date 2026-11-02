"""Five-display-item manuscript proposal from verified aggregate results only."""
import os,csv,json,hashlib,argparse,shutil,subprocess,sys
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/clif_manuscript_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,PathPatch,Rectangle
from matplotlib.path import Path as MPath
from matplotlib.ticker import PercentFormatter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'output/run_manuscript_proposal_20261030');a=p.parse_args();OUT=a.output.resolve();OUT.mkdir(parents=True,exist_ok=False)
active=json.loads((ROOT/'config/active_analysis.json').read_text());DESC=ROOT/active['exact_root']/'reports/publication_tables';COX=ROOT/active['four_endpoint_cox_reports'];FED=ROOT/active['primary_root'];SANK=ROOT/active['sankey_root'];POOL=ROOT/active['pooled_survival_root']/'reports';ROC=ROOT/active['pooled_roc_root']/'reports'
inputs=set()
def read(p):
 inputs.add(p)
 with p.open() as f:return list(csv.DictReader(f))
def verified(p):
 inputs.add(p);m=json.loads(p.read_text());assert m['status']=='passed';return m
consensus=active['primary_consensus_id']
for folder in [COX,POOL,ROC]:
 m=verified(folder/('manifest.json' if folder!=POOL else 'pooling_manifest.json'))
 assert m.get('primary_consensus_id',m.get('consensus_id'))==consensus
assert verified(SANK/'manifest.json')['consensus_id']==consensus
cont=read(DESC/'pooled_continuous.csv');cat=read(DESC/'pooled_categorical.csv');den=read(DESC/'denominators_by_site.csv')
GROUPS=['Overall','P1','P2','P3','P4','I'];N={g:sum(int(r['n']) for r in den if r['cohort']=='primary_phenotypes' and r['group']==g) for g in GROUPS};assert N['Overall']==sum(N[g] for g in GROUPS[1:])==3310
C={(r['group'],r['variable']):r for r in cont if r['cohort']=='primary_phenotypes'}
K={}
for r in cat:
 if r['cohort']=='primary_phenotypes':K.setdefault((r['group'],r['variable']),[]).append(r)
def continuous(g,v):
 r=C[g,v];assert r['status'] in ['observed_all_sites','all_missing']
 if not int(r['observed_n']):return 'Unavailable'
 sd=f'{float(r["sd"]):.1f}' if r['sd'] else 'NA'
 return f'{float(r["mean"]):.1f} ({sd})'
def pct(n,g):return f'{n:,} ({100*n/N[g]:.1f}%)'
def category(g,v,values=None,nonfull=False):
 rows=K[g,v];assert all(r['status']=='observed_all_sites' for r in rows)
 assert sum(int(r['n']) for r in rows)==N[g]
 n=sum(int(r['n']) for r in rows if (r['category'] not in ['full','unknown'] if nonfull else r['category'] in values))
 return pct(n,g)
entries=[]
def row(label,fn):entries.append([label]+[fn(g) for g in GROUPS])
row('Hospitalizations, N',lambda g:f'{N[g]:,}')
row('Age, years',lambda g:continuous(g,'age_years'))
row('Female sex',lambda g:category(g,'sex_group',['Female']))
row('Black race',lambda g:category(g,'race_group',['Black']))
row('Race unknown',lambda g:category(g,'race_group',['Unknown']))
row('Noncancer Charlson index',lambda g:continuous(g,'charlson_noncancer'))
row('Advanced-cancer code proxy',lambda g:category(g,'advanced_cancer',['1','TRUE']))
row('Admission code status: full',lambda g:category(g,'code_at_admission',['full']))
row('Admission code status: nonfull',lambda g:category(g,'code_at_admission',nonfull=True))
row('Admission code status: unknown',lambda g:category(g,'code_at_admission',['unknown']))
row('Modified SOFA available',lambda g:pct(int(C[g,'sofa_total']['observed_n']),g))
row('Modified SOFA, observed cases',lambda g:continuous(g,'sofa_total'))
row('Adapted APACHE available',lambda g:pct(int(C[g,'apache_adapted_total']['observed_n']),g))
row('Adapted APACHE, observed cases',lambda g:continuous(g,'apache_adapted_total'))
row('Documented IMV during first 72 h',lambda g:category(g,'any_imv',['1','TRUE']))
row('Documented HFNC during first 72 h',lambda g:category(g,'any_hfnc',['1','TRUE']))
row('Documented vasopressor, first 72 h',lambda g:category(g,'any_vasopressor',['1','TRUE']))
row('Any positive hypoxemia evidence',lambda g:category(g,'any_hypox',['1','TRUE']))
row('Any positive hypercapnia evidence',lambda g:category(g,'any_hyper',['1','TRUE']))
row('Unknown-physiology hours / 72 h',lambda g:continuous(g,'physiology_unknown_hours'))
headers=['Characteristic','Overall','P1\nRoom air','P2\nLow-flow O2','P3\nHFNC','P4\nIMV','I\nInsufficient info']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
def table_png(name,title,subtitle,headers,body,widths,note,height):
 fig,ax=plt.subplots(figsize=(13,height));ax.axis('off');fig.subplots_adjust(left=.035,right=.98,top=.82,bottom=.20)
 fig.suptitle(title,x=.035,ha='left',fontsize=20,y=.97);fig.text(.035,.89,subtitle,fontsize=11)
 tab=ax.table(cellText=body,colLabels=headers,colWidths=widths,cellLoc='center',loc='center');tab.auto_set_font_size(False);tab.set_fontsize(9.2);tab.scale(1,1.62)
 for (i,j),cell in tab.get_celld().items():
  cell.set_edgecolor('#DCE4E9');cell.set_linewidth(.45)
  if i==0:cell.set_facecolor('#E5EDF3');cell.set_text_props(weight='bold');cell.set_height(cell.get_height()*1.7)
  else:cell.set_facecolor('#F6F8FA' if i%2==0 else 'white')
  if j==0:cell.set_text_props(ha='left');cell.PAD=.04
 fig.text(.035,.035,note,fontsize=9,color='#444444',va='bottom')
 fig.savefig(OUT/(name+'.png'),dpi=200,facecolor='white');plt.close(fig)
def csvwrite(name,headers,body):
 with (OUT/name).open('w') as f:w=csv.writer(f);w.writerow(headers);w.writerows(body)
table_png('table_1','Table 1. Characteristics of the landmark cohort','Pooled UCMC, NU and Rush; current three-site results, pending final site returns',headers,entries,[.27]+[.73/6]*6,
          'Values are mean (SD) or N (% of the column). Scores summarize observed cases; availability is shown separately. NA: SD not estimable for N=1.\n'
          'P1-P4 are shared support-ARF phenotypes; I indicates insufficient support information. Unknown status is retained.\n'
          'Advanced-cancer proxy: metastatic cancer, malignant pleural effusion or cachexia codes; it is not measured cancer stage.\n'
          'SOFA and APACHE are CLIF-derived adaptations, not validated standard scores. Site-specific characteristics and detailed missingness belong in the supplement.\n'
          'Support and vasopressor rows indicate documented use; zero does not establish absence when source coverage is incomplete.',10)
csvwrite('table_1.csv',headers,entries)
cox=read(COX/'table2_presentation.csv');body=[]
labels={'death':'Hospital death','hospice':'Hospice discharge','death_or_hospice':'Death or hospice','live_discharge':'Live non-hospice discharge*'}
for r in cox:
 body.append([labels[r['endpoint']],r['phenotype'],f'{int(r["cohort_events"]):,}',f'{float(r["hazard_ratio"]):.2f} ({float(r["conf_low"]):.2f}-{float(r["conf_high"]):.2f})',f'{float(r["i2_percent"]):.1f}%',r['contributing_sites']])
table2_headers=['Endpoint','Phenotype\nvs P1','Total\nevents','Adjusted HR (95% CI)','I-squared','Sites']
table_png('table_2','Table 2. Post-landmark outcome associations','Cause-specific Cox; 3,310 hospitalizations; follow-up to 28 days after ICU +72 hours',table2_headers,body,[.28,.12,.12,.26,.12,.10],
          'Core adjustment: age, sex, race, admission year, noncancer Charlson, advanced-cancer proxy and admission code status.\n'
          'Reference: P1. Patient-clustered SEs within site; REML meta-analysis with modified Hartung-Knapp 95% CIs.\n'
          'Competing hospital exits are censored. I is retained with a separate baseline hazard; no clinical HR is reported for I. Event totals repeat across contrasts.\n'
          '*Live discharge excludes hospice. P4 discharge hazards vary over time at all three sites; retain exploratory 0-7/7-28-day results in the supplement.',8)
csvwrite('table_2.csv',table2_headers,body)
# Figure 1: support and ARF together; Sankey standalone.
for source in [ROOT/active['exact_root']/'reports/aggregate_collection/exact_counts/profiles/hourly_support.csv',ROOT/active['exact_root']/'reports/aggregate_collection/exact_counts/profiles/hourly_joint.csv',SANK/'pooled_initial_support_to_phenotype_counts.csv']:inputs.add(source)
subprocess.run([sys.executable,str(ROOT/'diagnostics/arf_display/main_figure.py'),str(OUT/'figure_1.png')],check=True)
# Redraw verified aggregate curves without titles, subtitles or captions; no refitting.
survival_root=OUT/'survival';survival_reports=survival_root/'reports';survival_reports.mkdir(parents=True)
for name in ['pooled_survival.csv','site_survival.csv','pooled_discharge_incidence.csv','site_discharge_incidence.csv','pooling_manifest.json']:
 source=POOL/name;inputs.add(source);shutil.copy2(source,survival_reports/name)
subprocess.run([sys.executable,str(ROOT/'diagnostics/pooled_survival/plot.py')],env={**os.environ,'CLIF_POOLED_SURVIVAL_OUTPUT':str(survival_root)},check=True)
roc_reports=OUT/'roc'
shares=[ROOT/active['pooled_roc_root']/'sites'/site/'share' for site in ['UCMC','NU','RUSH']]
for share in shares:inputs.add(share/'manifest.json')
subprocess.run([sys.executable,str(ROOT/'diagnostics/pooled_roc/central_plot.py'),str(roc_reports),*[str(share) for share in shares]],check=True)
for source,target in [(survival_reports/'pooled_phenotype_survival.png','figure_2.png'),(roc_reports/'pooled_roc_death_or_hospice.png','figure_3.png')]:shutil.copy2(source,OUT/target)
for name in ['pooled_auc.csv','pooled_roc.csv','paired_auc_differences.csv']:
 assert (ROC/name).read_bytes()==(roc_reports/name).read_bytes(),f'ROC numerical output changed: {name}'
items=[('table_1','Table 1. Landmark cohort characteristics','Clinical and respiratory descriptors by shared phenotype. Detailed site tables and missingness will accompany the supplement.'),('table_2','Table 2. Four cause-specific outcomes','Death, hospice, death/hospice and live non-hospice discharge; core-adjusted pooled HRs. The discharge time-varying sensitivity remains essential.'),('figure_1','Figure 1. Respiratory support and ARF trajectories','Smoothed observed support and four-category ARF evidence in P1-P4; the Sankey is exported separately. The 130 insufficient-information stays remain in Table1 and the models.'),('figure_2','Figure 2. Post-landmark survival and discharge','Hospital-death cause-specific KM and Aalen-Johansen live non-hospice discharge incidence. Neither panel measures all-cause survival after discharge.'),('figure_3','Figure 3. Pooled held-out ROC at days 7, 14 and 28','Death/hospice. All four models use the same 341 complete-score stays: UCMC96, NU198, Rush47. This selected cohort differs from the 3,310-stay primary cohort.')]
items.append(('standalone_initial_support_sankey','Standalone Sankey (supplement proposal)','Initial respiratory support to shared phenotype; separate from the three main figures.'))
# Five main display items plus standalone Sankey, one per page.
pdf=ROOT/'output/pdf'/f'{OUT.name}_review.pdf';pdf.parent.mkdir(parents=True,exist_ok=True);c=canvas.Canvas(str(pdf),pagesize=(792,720));c.setTitle('CLIF lung cancer - proposed manuscript figures and tables')
for index,(stem,title,note) in enumerate(items,1):
 c.setFillColorRGB(.15,.22,.28);c.setFont('Helvetica-Bold',15);c.drawString(28,686,title)
 # Explanatory captions remain outside figure images.
 image=ImageReader(str(OUT/(stem+'.png')));iw,ih=image.getSize();w,h=736,570;scale=min(w/iw,h/ih);dw,dh=iw*scale,ih*scale;c.drawImage(image,28+(w-dw)/2,85+(h-dh)/2,width=dw,height=dh)
 c.setFillColorRGB(.25,.25,.25);c.setFont('Helvetica',9)
 import textwrap
 for line,j in zip(textwrap.wrap(note,125),range(3)):
  if stem.startswith('table_'):c.drawString(28,59-12*j,line)
 c.setFont('Helvetica',8);c.drawRightString(764,20,f'{index} / {len(items)}');c.showPage()
c.save()
html=['<!doctype html><meta charset="utf-8"><title>Proposed manuscript display items</title><style>body{font-family:Arial,sans-serif;max-width:1250px;margin:35px auto;color:#173044}img{width:100%;border:1px solid #ddd}section{margin-bottom:60px}p{line-height:1.5}</style><h1>Proposed manuscript: two tables and three figures</h1><p>Current three-site rehearsal. Final site returns and respiratory ETL corrections are pending. Primary N=3,310; matched-score ROC N=341.</p>']
for stem,title,note in items:html.append(f'<section><h2>{title}</h2><p>{note}</p><img src="{stem}.png"></section>')
(OUT/'index.html').write_text('\n'.join(html))
manifest=dict(status='passed',figure_captions_and_subtitles=False,standalone_sankey_in_main_item_count=False,primary_consensus_id=consensus,main_display_items=5,primary_n=N['Overall'],roc_n=341,patient_level_inputs_read=False,models_refitted=False,mean_SD='pooled exact sufficient statistics; no average site medians',counts='exact aggregate categorical partitions and denominators',table2='core cause-specific Cox, 28 days after landmark; discharge excludes hospice',figure1_smoothing='five-hour moving average of observed hourly proportions and convex display interpolation; definitions and assignments unchanged',figure1_shown_n=sum(N[g] for g in ['P1','P2','P3','P4']),figure1_omitted_from_display_only=N['I'],pending='Final site returns and respiratory ETL corrections',input_md5={str(p.relative_to(ROOT)):hashlib.md5(p.read_bytes()).hexdigest() for p in sorted(inputs)},files={p.name:hashlib.md5(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file()})
manifest['review_pdf']=str(pdf.relative_to(ROOT));manifest['review_pdf_md5']=hashlib.md5(pdf.read_bytes()).hexdigest()
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2));print(OUT)
