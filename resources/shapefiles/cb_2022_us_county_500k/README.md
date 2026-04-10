Bundled U.S. Census Bureau county shapefile used for site-level CONUS maps.

Source:
- https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_us_county_500k.zip

Notes:
- This folder is used by `/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/02_federated_clusters.R` before falling back to `tigris::counties()`.
- Keeping the shapefile in the repo allows sites to generate maps without requiring internet access.
