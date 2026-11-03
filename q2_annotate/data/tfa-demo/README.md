# TFA example taxonomy

The resolved lineages use these seven standard ranks:
domain (`d`), phylum (`p`), class (`c`), order (`o`), family (`f`), genus (`g`),
and species (`s`). NCBI's additional kingdoms, unranked clades, species complexes,
and strains are omitted so each numbered selector level has the same meaning.
Names and lineages were checked against NCBI Taxonomy on 2026-10-14.

| Example ID | Example assignment | NCBI lineage reference |
| --- | --- | --- |
| T001 | Escherichia coli | [ncbitaxon:562](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=562&mode=Info) |
| T002 | Klebsiella pneumoniae | [ncbitaxon:573](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=573&mode=Info) |
| T003 | Enterococcus faecium | [ncbitaxon:1352](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1352&mode=Info) |
| T004 | Bacteroides fragilis | [ncbitaxon:817](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=817&mode=Info) |
| T005 | Phocaeicola vulgatus | [ncbitaxon:821](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=821&mode=Info) |
| T006 | Faecalibacterium prausnitzii | [ncbitaxon:853](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=853&mode=Info) |
| T007 | Roseburia intestinalis | [ncbitaxon:166486](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=166486&mode=Info) |
| T008 | Akkermansia muciniphila | [ncbitaxon:239935](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=239935&mode=Info) |
| T009 | Bifidobacterium longum | [ncbitaxon:216816](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=216816&mode=Info) |
| T010 | Segatella (genus only) | [ncbitaxon:165179](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=165179&mode=Info) |
| T011 | Clostridioides difficile | [ncbitaxon:1496](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1496&mode=Info) |
| T012 | Rikenellaceae (family only) | [ncbitaxon:28117](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=28117&mode=Info) |
| T013 | Unclassified | — |

`T005` uses *Phocaeicola vulgatus* (formerly *Bacteroides vulgatus*).
Ten taxa are resolved to species; `T010` is deliberately resolved only to genus
(*Segatella*) and `T012` only to family (Rikenellaceae). Their linked species
records document the known ancestral lineage, not a species assignment in the
example. `T013` is explicitly `d__Unclassified`. These cases demonstrate rank
fallback labels and the Show unclassified taxa control. IDs, sample metadata,
gene inventory, and abundance values are unchanged.
