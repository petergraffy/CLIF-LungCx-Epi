# CLIF Consortium Email Builder

This Shiny app reads the consortium member spreadsheet and helps you:

- filter members by site
- narrow by member status
- search by name, email, or affiliation
- choose exactly which people to include
- copy a semicolon-separated or comma-separated recipient list
- export the current selection as CSV

## Run it

From the project root:

```r
shiny::runApp("tools/clif_member_email_app")
```

The app defaults to:

`/Users/saborpete/Downloads/Untitled spreadsheet.xlsx`

You can also upload a different workbook from inside the app.

## Notes

- Site labels are inferred from affiliation text and email domains.
- If `Gmail` is available for a member, you can choose to prefer that over the work email.
- Missing `Member status` values are shown as `Unknown`.
