Four post-landmark cause-specific Cox endpoints

Endpoints:
- death: hospital death; hospice and live discharge are competing exits.
- hospice: hospice discharge; hospital death and other live discharge compete.
- death_or_hospice: composite; other live discharge competes.
- live_discharge: live NON-HOSPICE discharge; death and hospice compete.
Competing exits are censored at their observed exit time in cause-specific Cox.
Follow-up is capped at28days after each analysis-specific ICU landmark. A live-
discharge HR >1 means a higher instantaneous discharge rate among stays still
hospitalized; it is not an absolute probability or the inverse of a death HR.
Retain Aalen-Johansen cumulative incidence for absolute discharge probabilities.

Future fresh primary/sensitivity outcome passes automatically fit all four via
fed_cox_associations() in diagnostics/federated/cox_lib.R, called by both workers.
Unadjusted, core and measured-SOFA sensitivity specifications use identical
endpoint cohorts within a specification. Core: age, sex, race, admission year,
noncancer Charlson, advanced-cancer proxy and initial code status; P1 reference.
All stays retained; I has a separate baseline hazard and no clinical contrast.
Patient-clustered standard errors and Efron ties. Warnings/nonfinite coefficients
are marked nonestimable, with no fabricated zero effects. Endpoint-specific
covariance, N/event/status and PH diagnostic exports accompany associations.
The generic earlier raw-profile workflow also fits all four but preserves its
existing full-hospital follow-up specification; it is not the adopted primary.

Completed frozen analyses need no rediscovery or prediction refit. Run locally:
Rscript diagnostics/four_endpoint_cox/site_capture.R SITE PRIMARY REVISION NEW_OUTPUT
This refits four endpoints from saved full-version assignments for primary and
all8 inferential sensitivities. Original3 endpoint coefficients and clustered
covariances must reproduce, assignments/cache hashes must remain unchanged.
PH reconstruction and exploratory0-7/7-28day effects include all four endpoints.
Return NEW_OUTPUT/share/ only; models/data remain under local/.

Coordinator:
Rscript diagnostics/four_endpoint_cox/collect.R NEW_REPORT SITE_SHARE [SITE_SHARE ...]
Only aggregate returns are read. Site coefficients/SEs are pooled by REML and
modified Hartung-Knapp, requiring >=2 estimable sites. Unavailable effects are
reported explicitly. Both overall and exploratory interval effects are pooled.
No new outcome predictors or ROC markers are trained by this extension.

Append to a complete final site packet:
Rscript diagnostics/four_endpoint_cox/add_to_return.R SITE FOUR_COX_SHARE OLD_RETURN NEW_RETURN
Exports live under federated/four_endpoint_cox/, with four_endpoint_cox_complete
and an explicit primary-export path in the site manifest. Earlier primary files
remain historical; use this extension for the updated four-endpoint Table2.
CLIF_REQUIRE_FOUR_ENDPOINT_COX=true makes the site_return collector require the
extension at every site. Mixed endpoint phases/consensus IDs are rejected.

Three-site rehearsal: bash diagnostics/four_endpoint_cox/run.sh
Set CLIF_FOUR_COX_PRIMARY, CLIF_FOUR_COX_REVISION and a fresh CLIF_FOUR_COX_OUTPUT.
Current completed extension: output/run_nearest_assignment_20261030/four_endpoint_cox_v2/
Current final report folder: reports_with_piecewise/
27 site-variant passes, 324 Cox fits; all previously fitted endpoint effects
reproduced. Table2: table2_four_endpoint_core_cox.csv in the final report folder.
Primary discharge PH checks flag P4 at all3sites; interpretation must retain the
exploratory interval results. These intervals were added after inspecting PH
checks and are not prespecified confirmatory time-varying effects.

Test: Rscript tests/four_endpoint_cox.R
Method reference: https://stat.ethz.ch/CRAN/web/packages/survivalVignettes/vignettes/tutorial.html
