# How to publish this repository and get the DOI

You need to do this yourself (it requires your GitHub and Zenodo accounts). The whole thing takes
about fifteen minutes. Zenodo is what gives you a permanent DOI, which is what the manuscript cites.

## 1. Replace the placeholder

`sessionInfo.txt` is a placeholder. In the R session that produced the final results:

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

Overwrite the file with that output. If you use renv, commit `renv.lock` as well. Do not publish with
the placeholder still in place - the manuscript tells readers that exact package versions are here.

## 2. Create the GitHub repository

- Suggested name: `appendiceal-nodal-risk-surgical-extent`
- Public, no template files, no README (this folder already has one)
- Upload the contents of this folder (README.md, DATA_DICTIONARY.md, LICENSE, CITATION.cff,
  .gitignore, sessionInfo.txt, `R/`, `results_reference/`)
- Check after uploading that no CSV of patient-level data went up; `.gitignore` blocks the usual
  names, but the check is worth thirty seconds

## 3. Mint the DOI on Zenodo

- Sign in at <https://zenodo.org> with the GitHub option
- Go to the GitHub tab, find the repository, and set its switch to **On**
- Back on GitHub: Releases, then "Create a new release", tag `v1.0.0`, title
  "Analysis code for the appendiceal adenocarcinoma nodal-risk study", publish
- Zenodo archives the release within a minute or two and issues a DOI. Use the **Concept DOI**
  (the one Zenodo labels "all versions") so that later releases stay covered by the same citation

## 4. Put the DOI into the submission package

Two places need the real identifier:

1. `01_Manuscript/Manuscript.docx`, Data Availability Statement - replace
   `[INSERT REPOSITORY URL/DOI]` with, for example:
   `https://doi.org/10.5281/zenodo.XXXXXXX (GitHub: https://github.com/<account>/appendiceal-nodal-risk-surgical-extent)`
2. `CITATION.cff` in the repository - replace `<your-account>` in `repository-code`

The cover letter already says the code is publicly available, so it needs no change once the DOI exists.

## If you would rather not use Zenodo

A plain public GitHub URL satisfies the journal. A DOI is better because it is permanent and citable,
and because MDPI's data-availability guidance asks for stable identifiers. If you use GitHub only,
cite the URL together with the commit hash of the released state.
