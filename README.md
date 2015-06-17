# aCGH_viewer

### To visualize aCGH profiles using a segmentation table of the form of the DNAcopy (CBS) output
### Files should look like this:

| ID | chrom | loc.start | loc.end | num.mark | seg.mean |
|----|-------|-----------|---------|----------|----------|
| Sample.1 | 1 | 882803 | 61849262 | 3475 | 0.0688 |
| Sample.1 | 1 | 61851376 | 61861128 | 2 | 1.3717 |
| Sample.1 | 1 | 61865995 | 249198060 | 11436 | 0.0013 |
| Sample.1 | 2 | 249266324 | 492024204 | 16500 | -0.1146 |
...

### Pay attention you run the segmentation step using the chromosomal probe locations.

### Two screens for a full information
1. The genomic profile
2. The gene values

### aCGH Viewer features:
1. Load a segmentation file stored locally
2. display genes using valid HUGO symbols
3. Select the full profil, or a specific chromosome
4. Change the gain/loss color palette
5. Recenter the entire profile
6. Rescale the y-axis
7. Define the Gain/Loss and segment length thresholds to only display the relevant regions
8. Download the profile as it appears on the screen
9. Download the gene values

### The gene values are automatically updated with any change on the profile.
