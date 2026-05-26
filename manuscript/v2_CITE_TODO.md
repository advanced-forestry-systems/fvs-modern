# v2 manuscript: unresolved `[CITE: ...]` markers

*Generated 2026-05-16 from a grep across the four v2 manuscript files.
Each row gives the marker location, what is being claimed, and the
suggested citation form. Resolves to a Paperpile / BibTeX entry list at
the bottom.*

## Counts

| File | Markers |
|---|---|
| `v2_abstract.md` | 3 (one is meta-commentary, not a real citation) |
| `v2_section1_introduction.md` | 6 |
| `v2_section5_discussion.md` | 6 (three are checklist items, not real citations) |
| `v2_section6_conclusion.md` | 0 |
| **Total real citations needed** | **9 distinct sources** |

## Markers in document order

### v2_abstract.md

| Line | Marker | Claim being supported | Resolution |
|---|---|---|---|
| 5 | `[CITE: ...]` (meta) | Documentation note about there being two markers | meta-text, not a real citation; can be removed in final |
| 50 | `[CITE: 10.5281/zenodo.19802673]` | Zenodo DOI for the fvs-modern release | DOI is correct; this should be a standard preprint-style citation: "fvs-modern v2026.05.x [Software], DOI: 10.5281/zenodo.19802673" |
| 59 | `[CITE: pre-specified threshold reference]` | The 5% functional-equivalence threshold | Needs a methodological citation. Suggest: pre-specified analysis plan / SAP, or Weiskittel et al. 2011 (if the threshold was first articulated there), or a methods paper on Bayesian equivalence testing. If no prior reference, add explicit sentence: "We defined the 5% threshold a priori before fitting B1, drawing on the equivalence-testing framework of [Lakens et al. 2018; Westlake 1981]." |

### v2_section1_introduction.md

| Line | Marker | Claim being supported | Resolution |
|---|---|---|---|
| 14 | `[CITE: Dixon 2002]` | FVS history and operational status | **Dixon, G.E. 2002.** Essential FVS: A User's Guide to the Forest Vegetation Simulator. USDA Forest Service, Internal Report. Fort Collins, CO. |
| 19 | `[CITE: Weiskittel et al. 2011]` | Pre-2005 inventory data underlying default parameters | **Weiskittel, A.R., Hann, D.W., Kershaw, J.A., Vanclay, J.K. 2011.** Forest Growth and Yield Modeling. Wiley-Blackwell. |
| 40 | `[CITE: Bechtold and Patterson 2005]` | FIA annual panel design | **Bechtold, W.A., Patterson, P.L. (Eds.) 2005.** The Enhanced Forest Inventory and Analysis Program — National Sampling Design and Estimation Procedures. USDA Forest Service GTR-SRS-80. |
| 41 | `[CITE: Gelman et al. 2013]` | Bayesian hierarchical modeling reference | **Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., Rubin, D.B. 2013.** Bayesian Data Analysis, 3rd ed. CRC Press. |
| 57 | `[CITE: Clark et al. 2016; Bohn and Huth 2018]` | Bayesian methods in tree growth modeling | **Clark, J.S. et al. 2016.** "Tree growth inference and prediction when forest data are heterogeneous." Forest Ecology and Management. **Bohn, F.J., Huth, A. 2018.** "The importance of forest structure to biodiversity-productivity relationships." Royal Society Open Science 5: 172395. |
| 71 | `[CITE: Woodall et al. 2025]` | Forest Carbon Modeling Group statement on uncertainty | **Woodall, C.W. et al. 2025.** Forest Carbon Modeling Group consensus statement. (Source TBD — likely a special-issue editorial in Forest Ecology and Management or Carbon Balance and Management.) |

### v2_section5_discussion.md

| Line | Marker | Claim being supported | Resolution |
|---|---|---|---|
| 8 | `[CITE: ...]` (meta) | Documentation note | meta-text, remove in final |
| 39 | `[CITE: pre-specified threshold reference]` | Same 5% threshold | Same resolution as Abstract line 59 |
| 125 | `[CITE: TRY]` | TRY plant trait database | **Kattge, J. et al. 2020.** "TRY plant trait database — enhanced coverage and open access." Global Change Biology 26: 119-188. https://doi.org/10.1111/gcb.14904 |
| 198 | `[CITE: figure numbers finalized at production]` | Internal cross-reference | Resolve during typesetting; not a real citation. Use "Figures 1, 3, 4" with whatever final numbering the typesetter uses. |
| 204 | `[CITE: pre-specified threshold reference]` (checklist) | Same threshold checklist item | tracked in resolution above |
| 205 | `[CITE: TRY]` (checklist) | Same TRY checklist item | tracked in resolution above |

## Consolidated reference list (9 sources)

Once the Paperpile entries are confirmed, replace each `[CITE: ...]`
marker with the in-text citation form preferred by the target journal
(e.g., Dixon 2002, Weiskittel et al. 2011, Bechtold and Patterson 2005,
Gelman et al. 2013, Clark et al. 2016, Bohn and Huth 2018,
Kattge et al. 2020, Woodall et al. 2025) and add the corresponding
BibTeX entries to `references.bib`. The Zenodo DOI is already an
artifact of the release pipeline; the in-text form is "fvs-modern
(v2026.05.x) [Software], doi:10.5281/zenodo.19802673".

The 5% functional-equivalence threshold needs a methods-paper anchor.
If none exists in the team's prior work, the cleanest resolution is to
add a single sentence near the threshold's first use: "We set the
threshold a priori at 5 percent based on similar equivalence-band
practice in [supporting reference]." Candidates for [supporting
reference] include the Lakens et al. 2018 TOST methods paper, or any
forest-growth-and-yield paper that has previously used the same
threshold.

## Suggested Paperpile pass workflow

1. Open each v2 .md file in turn (`v2_abstract.md`,
   `v2_section1_introduction.md`, `v2_section5_discussion.md`,
   `v2_section6_conclusion.md`).
2. For each `[CITE: ...]` marker, drag the corresponding Paperpile entry
   into the manuscript at the marker's location. The 9 sources above
   should all already be in your library; the cleanest is to drag them
   in via the Paperpile browser extension while editing the .docx, then
   propagate back to .md.
3. Replace the meta-text `[CITE: ...]` markers (Abstract line 5,
   Discussion line 8, Discussion line 198) with concrete text or
   delete the surrounding sentence.
4. Resolve the 5% threshold marker by adding a single new sentence or
   citing a prior team paper.
5. Run `pandoc fvs_combined_draft.md -o fvs_combined_draft.docx
   --reference-doc=fvs_combined_draft.docx` once the citations are
   in to regenerate the .docx for co-author distribution.
