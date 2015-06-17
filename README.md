# -> aCGH_viewer <-

#### Interactively visualize aCGH profiles from uploaded segmentation tables.

### ->Download an [example file](https://github.com/fredcommo/aCGH_viewer/tree/master/example), then [Try me](https://fredcommo.shinyapps.io/aCGH_viewer)<-

#### Inputs are segmentation tables (.csv or .tsv) of the same form as the [DNAcopy](http://www.bioconductor.org/packages/release/bioc/vignettes/DNAcopy/inst/doc/DNAcopy.pdf) (CBS) outputs.

#### Input format

| ID | chrom | loc.start | loc.end | num.mark | seg.mean |
|----|-------|-----------|---------|----------|----------|
| Sample.1 | 1 | 882803 | 118165973 | 30208 | 0.0283 |
| Sample.1 | 1 | 118166096 | 119534291 | 365 | 1.2952 |
| Sample.1 | 1 | 119535183 | 154853295 | 1868 | 2.5673 |
| ... | ... | ... | ... | ... | ... |
| Sample.1 | 2 | 15703 | 235214315 | 63763 | 0.8037 |
| Sample.1 | 2 | 235214807 | 242775910 | 2064 | 0.6753 |
| ... | ... | ... | ... | ... | ... |

#### Pay attention you run the segmentation step using the chromosomal probe locations.

#### Two screens for a full information
1. The genomic profile
2. The gene values

#### aCGH Viewer features:
1. Load a segmentation file stored locally
2. display genes using valid HUGO symbols
3. Select the full profil, or a specific chromosome
4. Change the gain/loss color palette
5. Recenter the entire profile
6. Rescale the y-axis
7. Define the Gain/Loss and segment length thresholds to only display the relevant regions
8. Download the profile as it appears on the screen
9. Download the gene values

#### The gene values are automatically updated with any change on the profile.

#### The confidentiality is preserved: uploaded files are instantaneously removed. Nothing is stored.

### Genomic plot
![alt tag](https://github.com/fredcommo/aCGH_viewer/blob/master/screenshots/screen1.png)

### Gene values
![alt tag](https://github.com/fredcommo/aCGH_viewer/blob/master/screenshots/screen2.png)

