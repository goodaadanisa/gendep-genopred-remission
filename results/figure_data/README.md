# Figure-source data

This directory contains **27 shareable source files** supporting the dissertation figures. Twenty-six files fully regenerate 12 figures. `supplementary_figure_s4_summary.tsv` provides a release-safe aggregate summary for Supplementary Figure S4; its participant-level ancestry-probability input remains controlled and is not included. The final S4 PNG and plotting code are retained.

## File allocation

| Figure | Source files | Purpose |
|---|---:|---|
| Figures 1, 2, 4 and 5 | 1 each | A single table contains all values required for each figure. |
| Figure 3 | 6 | Model summaries, ROC coordinates, repeat-level AUCs, AUC-increment intervals, calibration bins and formal calibration statistics use different data structures. |
| Supplementary Figure S1 | 2 | Variant-frequency plotting points and their aggregate summary. |
| Supplementary Figure S2 | 2 | Chromosome-level results and the genome-wide summary. |
| Supplementary Figure S3 | 3 | PRS correlations, distribution summaries and violin-density coordinates. |
| Supplementary Figure S4 | 1 | Public aggregate summary only; the participant-level plotting input remains controlled. |
| Supplementary Figures S5–S7 | 2 each | Plotting values accompanied by a concise provenance table. |
| Supplementary Figure S8 | 3 | Tuning summaries, outer-fit coefficient values and coefficient-stability summaries. |

Separate files are retained where figure components have different columns, units or levels of aggregation. This is clearer and more reusable than combining unrelated structures into one large table.

The committed source files contain no participant identifiers, individual genotype calls, clinical records, held-out participant predictions or participant-level PRS. `supplementary_figure_s1_frequency_points.tsv.gz` is a gzip-compressed three-column TSV containing 147,438 plotting rows. It can be read directly with standard R or Python tools without manual decompression.
