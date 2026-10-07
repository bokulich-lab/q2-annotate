# Synthetic TFA explorer dataset

This example is **simulated**, not derived from biological samples. It is
designed to exercise the TFA visualizer with 60 samples (30 paired subjects,
baseline and day 14), 13 taxon IDs, 14 functions, and 58 observed
taxon–function pairs. Fifteen subjects are assigned to each of two synthetic
treatment groups. The generator introduces a larger day 14 resistance signal
in selected taxa in the antibiotic group, along with lower simulated loads
for some commensal taxa. These patterns are illustration only and should not
be interpreted as biological findings.

`tfa-explorer-demo.qzv` is a ready-to-open visualization generated from the
artifacts and metadata below. Open it with `qiime tools view` or upload it to
QIIME 2 View.

For comparison, `tfa-explorer-demo-heatmap-metric.qzv` is a separately named
copy of the current visualization, and `tfa-explorer-demo-c47115a.qzv` was
generated with the same inputs from the earlier sum-only version.

The source files are readable TSVs:

- `tfa-table.tsv`: one row per taxon–function pair, one column per sample;
  zeros represent absent loads. The corresponding QIIME artifact is
  `feature-load.qza` (`FeatureTable[Frequency % Properties("tfa")]`).
- `gene-taxonomy.tsv`: maps each stable feature ID to its taxon ID, gene ID,
  and taxonomy label. The corresponding artifact is `gene-taxonomy.qza`
  (`FeatureData[Taxonomy % Properties('tfa')]`), required by the visualizer.
- `taxonomy.tsv`: source labels used by the generator to populate the
  gene-taxonomy mapping. `T013` has no simulated assignment and uses its ID
  as the mapping label.
- `sample-metadata.tsv`: subject, treatment, visit, combined treatment/visit,
  batch, and age. The categorical columns appear in the visualizer's group
  selector.

With this branch installed in a QIIME 2 environment, run from the repository
root:

```bash
qiime annotate explore-tfa \
  --i-feature-load examples/tfa_explorer/feature-load.qza \
  --i-gene-taxonomy examples/tfa_explorer/gene-taxonomy.qza \
  --m-metadata-file examples/tfa_explorer/sample-metadata.tsv \
  --o-visualization tfa-explorer-demo.qzv
```

Open `tfa-explorer-demo.qzv` with `qiime tools view`. Select
`treatment_visit` to compare the four groups in the boxplots. Choose Level 2
to combine taxa by phylum or Level 3 by genus. The taxon dropdown, heatmap,
and sample plots update together; Original taxa restores the individual taxa. Switch the
heatmap metric between sum, mean, and median to change its colors and link
ranking across all samples. Try `blaCTX-M` and
`butyryl_CoA_transferase` in the function selector, then select individual
taxa or click heatmap cells. Taxonomy labels come from the gene-taxonomy
mapping. `T013` retains its ID label from the mapping.

To regenerate all TSV, QZA, and current QZV files from the fixed random seed, run
`python examples/tfa_explorer/generate_demo.py` in the same environment.

Keep the frequency table and gene-taxonomy mapping together. The historical
`c47115a` visualization is self-contained and remains available for comparison;
its earlier input artifacts used the retired `TFA` type.
