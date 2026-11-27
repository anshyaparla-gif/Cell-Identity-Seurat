# GitHub Upload Checklist

Use this checklist before pushing the project to GitHub.

## What to commit

- [ ] [Scripts/](Scripts/) analysis scripts
- [ ] [figures/](figures/) final PNG and PDF figure outputs
- [ ] [results/](results/) small summary tables and CSV outputs
- [ ] [README.md](README.md)
- [ ] [project_context.md](project_context.md) if you want to preserve the working notes
- [ ] Any small helper files that are part of the workflow

## What to leave out

- [ ] Raw Xenium data archives and extracted raw data folders
- [ ] Large temporary files such as `.RDataTmp*`
- [ ] `.Rproj.user/`
- [ ] `.Rhistory` files if they contain messy console history you do not want public
- [ ] Large intermediate R objects that are not needed to rerun the analysis
- [ ] Scratch plots like `Rplot*.pdf` and `Rplot.png` unless they are final figures

## Recommended repository layout

- `Scripts/` = one clean script per figure or analysis
- `figures/` = final plot images and PDFs
- `results/` = summary tables and CSV outputs
- `objects/` = saved analysis objects that are small enough to share
- `Xenium_data/` = local raw data only, not tracked in GitHub

## Before the first push

- [ ] Run the scripts once to confirm they still work
- [ ] Check that filenames are clear and match the paper figures
- [ ] Verify that no large raw-data files are staged for commit
- [ ] Review `git status` before pushing
- [ ] Add a short README note describing the project and the figures

## Good rule of thumb

If a file is needed to **reproduce the final paper figures or explain the project**, keep it.
If a file is only an intermediate scratch artifact, leave it out.
