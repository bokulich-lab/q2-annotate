from importlib import resources

import pandas as pd
from biom import Table
from rachis import Artifact, Metadata
from scipy.sparse import csr_matrix

_EXAMPLE_DATA = resources.files("q2_annotate") / "data" / "tfa-demo"


def _table_from_example(filename):
    """Import the readable feature-by-sample TSV as a sparse artifact."""
    with (_EXAMPLE_DATA / filename).open() as source:
        frame = pd.read_csv(source, sep="\t", index_col=0)
    return Artifact.import_data(
        "FeatureTable[Frequency]",
        Table(csr_matrix(frame.to_numpy()), frame.index, frame.columns),
    )


def factory_contig_abundance():
    """Load the synthetic contig-by-sample abundance table."""
    return _table_from_example("contig-abundance.tsv")


def factory_feature_inventory():
    """Load the synthetic gene-by-contig count table."""
    return _table_from_example("gene-inventory.tsv")


def factory_taxonomy():
    """Load the example taxon labels, including unresolved assignments."""
    with (_EXAMPLE_DATA / "taxonomy.tsv").open() as source:
        taxonomy = pd.read_csv(source, sep="\t", index_col=0, dtype=str)
    return Artifact.import_data("FeatureData[Taxonomy]", taxonomy)


def factory_taxonomy_map():
    """Load the example mapping from taxon IDs to their contigs."""
    with (_EXAMPLE_DATA / "taxon-to-contigs.tsv").open() as source:
        assignments = pd.read_csv(source, sep="\t", dtype=str)
    mapping = assignments.groupby("Taxon ID")["Contig ID"].apply(list).to_dict()
    return Artifact.import_data("FeatureMap[TaxonomyToContigs]", mapping)


def factory_sample_metadata():
    """Load subject, treatment, visit, and other sample annotations."""
    with resources.as_file(_EXAMPLE_DATA / "sample-metadata.tsv") as path:
        return Metadata.load(path)


def estimate_tfa_example(use):
    """Explore simulated taxon/function links in a paired, 60-sample cohort."""
    use.comment(
        "Simulated data: 60 samples from 30 paired subjects, 13 taxon IDs, "
        "14 functions, and 58 observed taxon/function pairs."
    )
    contig_abundance = use.init_artifact("contig_abundance", factory_contig_abundance)
    feature_inventory = use.init_artifact(
        "feature_inventory", factory_feature_inventory
    )
    taxonomy = use.init_artifact("taxonomy", factory_taxonomy)
    taxon_to_contig_map = use.init_artifact("taxon_to_contig_map", factory_taxonomy_map)

    sample_metadata = use.init_metadata("sample_metadata", factory_sample_metadata)

    tfa, gene_taxonomy, visualization = use.action(
        use.UsageAction(plugin_id="annotate", action_id="estimate_tfa"),
        use.UsageInputs(
            contig_abundance=contig_abundance,
            feature_inventory=feature_inventory,
            taxonomy=taxonomy,
            taxon_to_contig_map=taxon_to_contig_map,
            metadata=sample_metadata,
        ),
        use.UsageOutputNames(
            tfa="tfa", gene_taxonomy="gene_taxonomy", visualization="visualization"
        ),
    )

    tfa.assert_output_type("FeatureTable[Frequency % Properties('tfa')]")
    gene_taxonomy.assert_output_type("FeatureData[Taxonomy % Properties('tfa')]")
    visualization.assert_output_type("Visualization")

    use.comment(
        "Open the visualization and group samples by treatment_visit to compare "
        "the four treatment/visit groups. These data are synthetic and do not "
        "represent biological findings."
    )
