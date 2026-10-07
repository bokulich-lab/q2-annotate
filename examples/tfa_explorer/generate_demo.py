"""Create a synthetic gut-microbiome TFA dataset and QIIME 2 artifacts."""

from pathlib import Path

import biom
import numpy as np
import pandas as pd
import qiime2
from qiime2.sdk import PluginManager
from scipy.sparse import csr_matrix

from q2_annotate.plugin_setup import plugin as annotate_plugin
from q2_annotate.tfa.tfa import _gene_taxonomy_id
from q2_types.plugin_setup import plugin as types_plugin

HERE = Path(__file__).resolve().parent
TAXA = {
    "T001": ("Pseudomonadota", "Escherichia", "coli"),
    "T002": ("Pseudomonadota", "Klebsiella", "pneumoniae"),
    "T003": ("Bacillota", "Enterococcus", "faecium"),
    "T004": ("Bacteroidota", "Bacteroides", "fragilis"),
    "T005": ("Bacteroidota", "Bacteroides", "vulgatus"),
    "T006": ("Bacillota", "Faecalibacterium", "prausnitzii"),
    "T007": ("Bacillota", "Roseburia", "intestinalis"),
    "T008": ("Verrucomicrobiota", "Akkermansia", "muciniphila"),
    "T009": ("Actinomycetota", "Bifidobacterium", "longum"),
    "T010": ("Bacteroidota", "Prevotella", "copri"),
    "T011": ("Bacillota", "Clostridioides", "difficile"),
    "T012": ("Bacteroidota", "Alistipes", "putredinis"),
    "T013": None,  # missing assignment exercises the ID fallback
}
PAIRS = {
    "T001": "carbohydrate_transport folate_biosynthesis blaTEM blaCTX-M qnrS sul1",
    "T002": "carbohydrate_transport folate_biosynthesis blaTEM blaCTX-M qnrS sul1",
    "T003": "carbohydrate_transport amino_acid_fermentation tetM ermB vanA",
    "T004": "beta_glucosidase mucin_glycosidase carbohydrate_transport folate_biosynthesis tetM",
    "T005": "beta_glucosidase mucin_glycosidase carbohydrate_transport tetM",
    "T006": "butyryl_CoA_transferase butyrate_kinase beta_glucosidase folate_biosynthesis",
    "T007": "butyryl_CoA_transferase butyrate_kinase beta_glucosidase carbohydrate_transport",
    "T008": "mucin_glycosidase beta_glucosidase folate_biosynthesis carbohydrate_transport",
    "T009": "beta_glucosidase carbohydrate_transport folate_biosynthesis tetM",
    "T010": "beta_glucosidase carbohydrate_transport folate_biosynthesis amino_acid_fermentation",
    "T011": "amino_acid_fermentation carbohydrate_transport tetM ermB vanA",
    "T012": "beta_glucosidase carbohydrate_transport folate_biosynthesis amino_acid_fermentation",
    "T013": "carbohydrate_transport blaTEM sul1",
}
RESISTANCE = {"blaTEM", "blaCTX-M", "tetM", "ermB", "vanA", "qnrS", "sul1"}
OPPORTUNISTS = {"T001", "T002", "T003", "T011", "T013"}


def write_sources():
    rng = np.random.default_rng(20260925)
    samples = []
    for subject in range(1, 31):
        treatment = "antibiotic" if subject <= 15 else "placebo"
        age = int(rng.integers(22, 72))
        subject_scale = float(rng.lognormal(0, 0.25))
        for visit in ("baseline", "day14"):
            samples.append(
                {
                    "#SampleID": f"S{subject:02d}_{'BL' if visit == 'baseline' else 'D14'}",
                    "subject": f"P{subject:02d}",
                    "treatment": treatment,
                    "visit": visit,
                    "treatment_visit": f"{treatment}_{visit}",
                    "batch": "B1" if subject % 2 else "B2",
                    "age": age,
                    "scale": subject_scale,
                }
            )

    metadata = pd.DataFrame(samples).drop(columns="scale")
    with (HERE / "sample-metadata.tsv").open("w", encoding="utf-8") as output:
        output.write("\t".join(metadata.columns) + "\n")
        output.write(
            "#q2:types\t" + "\t".join(["categorical"] * 5 + ["numeric"]) + "\n"
        )
        metadata.to_csv(output, sep="\t", index=False, header=False)

    taxonomy = [
        {
            "Feature ID": taxon_id,
            "Taxon": f"d__Bacteria; p__{phylum}; g__{genus}; s__{genus} {species}",
        }
        for taxon_id, assignment in TAXA.items()
        if assignment is not None
        for phylum, genus, species in [assignment]
    ]
    pd.DataFrame(taxonomy).to_csv(HERE / "taxonomy.tsv", sep="\t", index=False)

    rows = []
    for taxon_id, functions in PAIRS.items():
        for function_id in functions.split():
            row = {"taxon_id": taxon_id, "function_id": function_id}
            pair_scale = float(rng.lognormal(0, 0.35))
            for sample in samples:
                exposed = (
                    sample["treatment"] == "antibiotic" and sample["visit"] == "day14"
                )
                opportunist = taxon_id in OPPORTUNISTS
                prevalence = 0.21 if opportunist else 0.34
                scale = pair_scale * sample["scale"]
                if function_id in RESISTANCE:
                    prevalence *= 0.68
                    scale *= 0.72
                if exposed and opportunist:
                    prevalence *= 1.75
                    scale *= 2.5
                    if function_id in RESISTANCE:
                        scale *= 2.1
                if exposed and not opportunist:
                    prevalence *= 0.65
                    scale *= 0.48
                row[sample["#SampleID"]] = (
                    round(float(rng.lognormal(np.log(7 * scale), 0.55)), 2)
                    if rng.random() < min(prevalence, 0.94)
                    else 0
                )
            rows.append(row)
    pd.DataFrame(rows).to_csv(HERE / "tfa-table.tsv", sep="\t", index=False)


def write_artifacts():
    # Register the checkout so this also works before installing the plugin.
    manager = PluginManager(add_plugins=False)
    manager.add_plugin(types_plugin, package="q2_types", project_name="q2-types")
    manager.add_plugin(
        annotate_plugin, package="q2_annotate", project_name="q2-annotate"
    )

    frame = pd.read_csv(HERE / "tfa-table.tsv", sep="\t")
    sample_ids = list(frame.columns[2:])
    table = biom.Table(
        csr_matrix(frame[sample_ids].to_numpy(dtype=float)),
        observation_ids=[
            _gene_taxonomy_id(taxon, function)
            for taxon, function in zip(frame.taxon_id, frame.function_id)
        ],
        sample_ids=sample_ids,
    )
    feature_load = qiime2.Artifact.import_data(
        "FeatureTable[Frequency % Properties('tfa')]", table
    )
    feature_load.save(HERE / "feature-load.qza")
    taxonomy = pd.read_csv(HERE / "taxonomy.tsv", sep="\t", index_col=0)
    taxonomy.index.name = "Feature ID"
    mapping = pd.DataFrame(
        {
            "Taxon": [
                taxonomy.at[taxon, "Taxon"] if taxon in taxonomy.index else taxon
                for taxon in frame.taxon_id
            ],
            "Taxon ID": frame.taxon_id.to_numpy(),
            "Gene ID": frame.function_id.to_numpy(),
        },
        index=pd.Index(table.ids(axis="observation"), name="Feature ID"),
    )
    mapping.to_csv(HERE / "gene-taxonomy.tsv", sep="\t")
    gene_taxonomy = qiime2.Artifact.import_data(
        "FeatureData[Taxonomy % Properties('tfa')]", mapping
    )
    gene_taxonomy.save(HERE / "gene-taxonomy.qza")
    result = annotate_plugin.visualizers["explore_tfa"](
        feature_load=feature_load,
        gene_taxonomy=gene_taxonomy,
        metadata=qiime2.Metadata.load(HERE / "sample-metadata.tsv"),
    )
    result.visualization.save(HERE / "tfa-explorer-demo.qzv")
    result.visualization.save(HERE / "tfa-explorer-demo-heatmap-metric.qzv")


if __name__ == "__main__":
    write_sources()
    write_artifacts()
