Rush respiratory-support coverage figure

Rscript diagnostics/rush_respiratory_gap/audit.R
python3 diagnostics/rush_respiratory_gap/plot.py

Output: output/run_RUSH_respiratory_gap_20261030/share/
The PNG/PDF and supporting CSVs/text contain Rush aggregate counts only.
The audit reads only the Rush frozen primary and source caches used in the
current three-site analysis. It does not read a newer release or the source EHR.

Findings:969 landmark stays;126 have no known support hour after waterfall.
Among those126,125 have no respiratory rows in the first72h or anywhere in the
supplied hospitalization; one has94 rows, all without device category. All126
have SpO2 and laboratory data in the window.118 have>=24 closed ICU-ADT hours.
Median closed ICU exposure67.0h (IQR43.6–72.0);47 cover the entire72h.
Annual denominators are printed separately for each supplied admission year.

This is evidence of missing respiratory-table coverage, not proof of an ETL
error: check room-air/off-oxygen extraction, flowsheet selection, admission joins
and timing against source records. Do not label undocumented support as normal.
No patient identifiers or other sites' data are in the share bundle.
