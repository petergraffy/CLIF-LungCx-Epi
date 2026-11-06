# Stage 2: one run, one return

Follow [the complete site instructions](../../code/README.md). Restore the pinned
R/Python environments and configure the source CLIF table paths. No Stage 1 folders or caches are needed. From the
repository root:

```sh
Rscript --vanilla scripts/site.R stage2-preflight config/config.json
Rscript --vanilla scripts/site.R stage2 config/config.json
```

Return only the printed `stage2_return/` directory. No second site pass is needed.
