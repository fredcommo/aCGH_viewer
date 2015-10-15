# aCGH_viewer

#### Interactively visualize aCGH profiles from uploaded segmentation tables.

### Download an [example file](https://github.com/fredcommo/aCGH_viewer/tree/master/example), then [Try me](https://fredcommo.shinyapps.io/aCGH_viewer)

#### Inputs are segmentation tables (.csv or .tsv) of the same form as the [DNAcopy](http://www.bioconductor.org/packages/release/bioc/vignettes/DNAcopy/inst/doc/DNAcopy.pdf) (CBS) outputs.

#### Minimal input format

| ID | chrom | loc.start | loc.end | num.mark | seg.mean |
|----|-------|-----------|---------|----------|----------|
| Sample.1 | 1 | 882803 | 118165973 | 30208 | 0.0283 |
| Sample.1 | 1 | 118166096 | 119534291 | 365 | 1.2952 |
| Sample.1 | 1 | 119535183 | 154853295 | 1868 | 2.5673 |
| ... | ... | ... | ... | ... | ... |
| Sample.1 | 2 | 15703 | 235214315 | 63763 | 0.8037 |
| Sample.1 | 2 | 235214807 | 242775910 | 2064 | 0.6753 |
| ... | ... | ... | ... | ... | ... |

#### Alternative input format, as provided by the [rCGH](http://bioconductor.org/packages/devel/bioc/html/rCGH.html) R package.

| ID | chrom | loc.start | loc.end | num.mark | seg.mean | probes.Sd |
|----|-------|-----------|---------|----------|----------|-----------|
| Sample.1 | 1 | 882803 | 118165973 | 30208 | 0.0283 | 0.8765 |
| Sample.1 | 1 | 118166096 | 119534291 | 365 | 1.2952 | 1.1234 |
| Sample.1 | 1 | 119535183 | 154853295 | 1868 | 2.5673 | 0.9765 |
| ... | ... | ... | ... | ... | ... | ... |
| Sample.1 | 2 | 15703 | 235214315 | 63763 | 0.8037 | 1.1348 |
| Sample.1 | 2 | 235214807 | 242775910 | 2064 | 0.6753 | 0.9583 |
| ... | ... | ... | ... | ... | ... | ... |

#### Important note:
##### Chromosomal probes locations must be used when the CNA object is built, not the genomic ones.

#### Two screens for a full information
1. The genomic profile
2. The gene values

#### aCGH Viewer features:
1. Load a segmentation file stored locally
2. Choose the genome build you are using
3. Display genes using valid HUGO symbols
4. Select the full profile, or isolate a specific chromosome
5. Change Gain/Loss colors
6. Merge segments shorter than a specified value, in Kb
7. Recenter the entire profile
8. Rescale the y-axis
9. Define the Gain/Loss and segment length thresholds to only display the relevant regions. The 'Genes table' will be filtered, accordingly.
10. Download the profile as it appears on the screen
11. Download the gene values

#### The gene values are automatically updated with any change made on the profile.

#### The confidentiality is preserved: uploaded files are instantaneously removed. Nothing is stored on the server.

### Genomic profile
![alt tag](https://github.com/fredcommo/aCGH_viewer/blob/master/screenshots/screen1.png)

### Gene values
![alt tag](https://github.com/fredcommo/aCGH_viewer/blob/master/screenshots/screen2.png)

#### Feel free to clone this repo, and to install the app on your own server.
