# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib
import platform

from q2_quality_control.plugin_setup import (
    filter_parameters,
    filter_parameter_descriptions,
)
from qiime2.plugin import Metadata
from q2_annotate.eggnog.types import (
    EggnogHmmerIdmapDirectoryFmt,
    EggnogHmmerIdmapFileFmt,
    EggnogHmmerIdmap,
)
from q2_types.bowtie2 import Bowtie2Index
from q2_types.profile_hmms import ProfileHMM, MultipleProtein, PressedProtein
from q2_types.feature_data import (
    FeatureData,
    Sequence,
    Taxonomy,
    ProteinSequence,
)
from q2_types.feature_table import (
    FeatureTable,
    Frequency,
    PresenceAbsence,
)
from q2_types.per_sample_sequences import (
    SequencesWithQuality,
    PairedEndSequencesWithQuality,
    JoinedSequencesWithQuality,
    MAGs,
    Contigs,
)
from q2_types.sample_data import SampleData
from q2_types.feature_map import FeatureMap, TaxonomyToContigs
from qiime2.core.type import (
    Bool,
    Range,
    Int,
    Str,
    Float,
    List,
    Choices,
    Visualization,
    Properties,
    TypeMap,
    TypeMatch,
    Threads,
)
from qiime2.plugin import Plugin, Citations
import q2_annotate._examples as ex
import q2_annotate
from q2_types.feature_data_mag import MAG
from q2_types.genome_data import NOG, Orthologs, GenomeData, Loci, Genes, Proteins
from q2_types.kaiju import KaijuDB
from q2_types.kraken2 import Kraken2Reports, Kraken2Outputs, Kraken2DB, Kraken2DBReport
from q2_types.kraken2 import BrackenDB
from q2_types.reference_db import (
    ReferenceDB,
    Diamond,
    Eggnog,
    NCBITaxonomy,
    EggnogProteinSequences,
)

citations = Citations.load("citations.bib", package="q2_annotate")

kraken2_params = {
    "threads": Int % Range(1, None),
    "confidence": Float % Range(0, 1, inclusive_end=True),
    "minimum_base_quality": Int % Range(0, None),
    "memory_mapping": Bool,
    "minimum_hit_groups": Int % Range(1, None),
    "quick": Bool,
    "report_minimizer_data": Bool,
}
kraken2_param_descriptions = {
    "threads": "Number of threads.",
    "confidence": "Confidence score threshold.",
    "minimum_base_quality": (
        "Minimum base quality used in classification. "
        "Only applies when reads are used as input."
    ),
    "memory_mapping": "Avoids loading the database into RAM.",
    "minimum_hit_groups": (
        "Minimum number of hit groups (overlapping "
        "k-mers sharing the same minimizer)."
    ),
    "quick": "Quick operation (use first hit or hits).",
    "report_minimizer_data": (
        "Include number of read-minimizers per-taxon and unique "
        "read-minimizers per-taxon in the report. If this parameter is "
        "enabled then merging kraken2 reports with the same sample ID from "
        "two or more input artifacts will not be possible."
    ),
}

partition_params = {"num_partitions": Int % Range(1, None)}
partition_param_descriptions = {
    "num_partitions": (
        "The number of partitions to split the contigs "
        "into. Defaults to partitioning into individual "
        "samples."
    )
}

plugin = Plugin(
    name="annotate",
    version=q2_annotate.__version__,
    website="https://github.com/bokulich-lab/q2-annotate",
    package="q2_annotate",
    description=(
        "MOdular SHotgun metagenome Pipelines with Integrated "
        "provenance Tracking: QIIME 2 plugin gor metagenome analysis with"
        "tools for genome binning and functional annotation."
    ),
    short_description="QIIME 2 plugin for metagenome analysis.",
)

importlib.import_module("q2_annotate.eggnog")

T_kraken_in_list, T_kraken_out_rep, T_kraken_out_hits = TypeMap(
    {
        List[
            SampleData[
                SequencesWithQuality
                | PairedEndSequencesWithQuality
                | JoinedSequencesWithQuality
            ]
        ]: (
            SampleData[Kraken2Reports % Properties("reads")],
            SampleData[Kraken2Outputs % Properties("reads")],
        ),
        List[SampleData[Contigs]]: (
            SampleData[Kraken2Reports % Properties("contigs")],
            SampleData[Kraken2Outputs % Properties("contigs")],
        ),
        List[FeatureData[MAG]]: (
            FeatureData[Kraken2Reports % Properties("mags")],
            FeatureData[Kraken2Outputs % Properties("mags")],
        ),
        List[SampleData[MAGs]]: (
            SampleData[Kraken2Reports % Properties("mags")],
            SampleData[Kraken2Outputs % Properties("mags")],
        ),
    }
)

plugin.pipelines.register_function(
    function=q2_annotate.kraken2.classification.classify_kraken2,
    inputs={
        "seqs": T_kraken_in_list,
        "db": Kraken2DB,
    },
    parameters={**kraken2_params, **partition_params},
    outputs=[
        ("reports", T_kraken_out_rep),
        ("outputs", T_kraken_out_hits),
    ],
    input_descriptions={
        "seqs": (
            "Sequences to be classified. Single-/paired-end reads,"
            "contigs, or assembled MAGs, can be provided."
        ),
        "db": "Kraken 2 database.",
    },
    parameter_descriptions={
        **kraken2_param_descriptions,
        **partition_param_descriptions,
    },
    output_descriptions={
        "reports": "Reports produced by Kraken2.",
        "outputs": "Output files produced by Kraken2.",
    },
    name="Perform taxonomic classification of reads or MAGs using Kraken 2.",
    description=(
        "Use Kraken 2 to classify provided DNA sequence reads, "
        "contigs, or MAGs into taxonomic groups."
    ),
    citations=[citations["wood2019"]],
)

T_kraken_in, T_kraken_out_rep, T_kraken_out_hits = TypeMap(
    {
        SampleData[
            SequencesWithQuality
            | PairedEndSequencesWithQuality
            | JoinedSequencesWithQuality
        ]: (
            SampleData[Kraken2Reports % Properties("reads")],
            SampleData[Kraken2Outputs % Properties("reads")],
        ),
        SampleData[Contigs]: (
            SampleData[Kraken2Reports % Properties("contigs")],
            SampleData[Kraken2Outputs % Properties("contigs")],
        ),
        FeatureData[MAG]: (
            FeatureData[Kraken2Reports % Properties("mags")],
            FeatureData[Kraken2Outputs % Properties("mags")],
        ),
        SampleData[MAGs]: (
            SampleData[Kraken2Reports % Properties("mags")],
            SampleData[Kraken2Outputs % Properties("mags")],
        ),
    }
)
plugin.methods.register_function(
    function=q2_annotate.kraken2.classification._classify_kraken2,
    inputs={
        "seqs": T_kraken_in,
        "db": Kraken2DB,
    },
    parameters=kraken2_params,
    outputs=[
        ("reports", T_kraken_out_rep),
        ("outputs", T_kraken_out_hits),
    ],
    input_descriptions={
        "seqs": (
            "Sequences to be classified. Single-/paired-end reads, "
            "contigs, or assembled MAGs can be provided."
        ),
        "db": "Kraken 2 database.",
    },
    parameter_descriptions=kraken2_param_descriptions,
    output_descriptions={
        "reports": "Reports produced by Kraken2.",
        "outputs": "Output files produced by Kraken2.",
    },
    name="Perform taxonomic classification of reads or MAGs using Kraken 2.",
    description=(
        "Use Kraken 2 to classify provided DNA sequence reads, "
        "contigs, or MAGs into taxonomic groups."
    ),
    citations=[citations["wood2019"]],
)

P_kraken_in, P_kraken_out = TypeMap(
    {
        Properties("reads", "contigs", "mags"): Properties("reads", "contigs", "mags"),
        Properties("reads", "contigs"): Properties("reads", "contigs"),
        Properties("reads", "mags"): Properties("reads", "mags"),
        Properties("contigs", "mags"): Properties("contigs", "mags"),
        Properties("reads"): Properties("reads"),
        Properties("contigs"): Properties("contigs"),
        Properties("mags"): Properties("mags"),
    }
)

plugin.methods.register_function(
    function=q2_annotate.kraken_helpers.collate_kraken2_reports,
    inputs={"reports": List[SampleData[Kraken2Reports % P_kraken_in]]},
    parameters={},
    outputs={"collated_reports": SampleData[Kraken2Reports % P_kraken_out]},
    name="Collate kraken2 reports.",
    description="Collates kraken2 reports.",
)

plugin.methods.register_function(
    function=q2_annotate.kraken_helpers.collate_kraken2_outputs,
    inputs={"outputs": List[SampleData[Kraken2Outputs % P_kraken_in]]},
    parameters={},
    outputs={"collated_outputs": SampleData[Kraken2Outputs % P_kraken_out]},
    name="Collate kraken2 outputs.",
    description="Collates kraken2 outputs.",
)

if platform.system() != "Darwin":
    plugin.methods.register_function(
        function=q2_annotate.kraken2.bracken.estimate_bracken,
        inputs={
            "kraken2_reports": SampleData[Kraken2Reports % Properties("reads")],
            "db": BrackenDB,
        },
        parameters={
            "threshold": Int % Range(0, None),
            "read_len": Int % Range(0, None),
            "level": Str % Choices(["D", "P", "C", "O", "F", "G", "S"]),
            "include_unclassified": Bool,
        },
        outputs=[
            ("reports", SampleData[Kraken2Reports % Properties("bracken")]),
            ("taxonomy", FeatureData[Taxonomy]),
            ("table", FeatureTable[Frequency]),
        ],
        input_descriptions={
            "kraken2_reports": "Reports produced by Kraken2.",
            "db": "Bracken database.",
        },
        parameter_descriptions={
            "threshold": (
                "Bracken: number of reads required PRIOR to abundance "
                "estimation to perform re-estimation."
            ),
            "read_len": (
                "Bracken: the ideal length of reads in your sample. "
                "For paired end data (e.g., 2x150) this should be set "
                "to the length of the single-end reads (e.g., 150)."
            ),
            "level": (
                "Bracken: specifies the taxonomic rank to  analyze. Each "
                "classification at this specified rank will receive an "
                "estimated number of reads belonging to that rank after "
                "abundance estimation."
            ),
            "include_unclassified": (
                "Bracken does not include the unclassified "
                "read counts in the feature table. Set this "
                "to True to include those regardless."
            ),
        },
        output_descriptions={
            "reports": "Reports modified by Bracken.",
        },
        name="Perform read abundance re-estimation using Bracken.",
        description=(
            "This method uses Bracken to re-estimate read abundances. "
            "Only available on Linux platforms."
        ),
        citations=[citations["wood2019"]],
    )

plugin.methods.register_function(
    function=q2_annotate.kraken2.build_kraken_db,
    inputs={"seqs": List[FeatureData[Sequence]]},
    parameters={
        "collection": Str
        % Choices(
            [
                "viral",
                "minusb",
                "standard",
                "standard8",
                "standard16",
                "pluspf",
                "pluspf8",
                "pluspf16",
                "pluspfp",
                "pluspfp8",
                "pluspfp16",
                "eupathdb",
                "nt",
                "corent",
                "gtdb",
                "greengenes",
                "rdp",
                "silva132",
                "silva138",
            ],
        ),
        "threads": Int % Range(1, None),
        "kmer_len": Int % Range(1, None),
        "minimizer_len": Int % Range(1, None),
        "minimizer_spaces": Int % Range(1, None),
        "no_masking": Bool,
        "max_db_size": Int % Range(0, None),
        "use_ftp": Bool,
        "load_factor": Float % Range(0, 1),
        "fast_build": Bool,
        "read_len": List[Int % Range(1, None)],
    },
    outputs=[
        ("kraken2_db", Kraken2DB),
        ("bracken_db", BrackenDB),
    ],
    input_descriptions={"seqs": "Sequences to be added to the Kraken 2 database."},
    parameter_descriptions={
        "collection": (
            "Name of the database collection to be fetched. "
            "Please check https://benlangmead.github.io/aws-"
            "indexes/k2 for the description of the available options."
        ),
        "threads": (
            "Number of threads. Only applicable when building a custom database."
        ),
        "kmer_len": "K-mer length in bp/aa.",
        "minimizer_len": "Minimizer length in bp/aa.",
        "minimizer_spaces": (
            "Number of characters in minimizer that are ignored in comparisons."
        ),
        "no_masking": (
            "Avoid masking low-complexity sequences prior to "
            "building; masking requires dustmasker or segmasker "
            "to be installed in PATH."
        ),
        "max_db_size": (
            "Maximum number of bytes for Kraken 2 hash table; "
            "if the estimator determines more would normally be "
            "needed, the reference library will be downsampled to fit."
        ),
        "use_ftp": "Use FTP for downloading instead of RSYNC.",
        "load_factor": "Proportion of the hash table to be populated.",
        "fast_build": (
            "Do not require database to be deterministically "
            "built when using multiple threads. This is faster, "
            "but does introduce variability in minimizer/LCA pairs."
        ),
        "read_len": (
            "Ideal read lengths to be used while building the Bracken database."
        ),
    },
    output_descriptions={
        "kraken2_db": "Kraken2 database.",
        "bracken_db": "Bracken database.",
    },
    name="Build Kraken 2 database.",
    description=(
        "This method builds Kraken 2 and Bracken databases either (1) "
        "from provided DNA sequences to build a custom database, or "
        "(2) simply fetches pre-built versions from an online "
        "resource."
    ),
    citations=[citations["wood2019"], citations["lu2017"]],
)

plugin.methods.register_function(
    function=q2_annotate.kraken2.database.inspect_kraken2_db,
    inputs={"db": Kraken2DB},
    parameters={"threads": Int % Range(1, None)},
    outputs=[("report", Kraken2DBReport)],
    input_descriptions={
        "db": "The Kraken 2 database for which to generate the report."
    },
    parameter_descriptions={"threads": "The number of threads to use."},
    output_descriptions={"report": "The report of the supplied database."},
    name="Inspect a Kraken 2 database.",
    description=(
        "This method generates a report of identical format to those "
        "generated by classify_kraken2, with a slightly different "
        "interpretation. Instead of reporting the number of inputs "
        "classified to a taxon/clade, the report displays the number "
        "of minimizers mapped to each taxon/clade."
    ),
    citations=[citations["wood2019"]],
)

select_features_taxonomy_description = (
    "Output taxonomy. Infra-clade ranks are ignored unless if they are "
    "strain-level. Missing internal ranks are annotated by their next "
    "most specific rank, with the exception of k__Bacteria and k__Archaea, "
    "which match their domain name."
)

select_features_description = (
    "Convert a Kraken 2 report, which is an annotated NCBI taxonomy tree, "
    "into generic artifacts for downstream analyses."
)

plugin.methods.register_function(
    function=q2_annotate.kraken2.kraken2_to_features,
    inputs={"reports": SampleData[Kraken2Reports]},
    parameters={"coverage_threshold": Float % Range(0, 100, inclusive_end=True)},
    outputs=[
        ("table", FeatureTable[PresenceAbsence]),
        ("taxonomy", FeatureData[Taxonomy]),
    ],
    input_descriptions={"reports": "Per-sample Kraken 2 reports."},
    parameter_descriptions={
        "coverage_threshold": (
            "The minimum percent coverage required to produce a feature."
        )
    },
    output_descriptions={
        "table": (
            "A presence/absence table of selected features. The features "
            "are not of even ranks, but will be the most specific rank "
            "available."
        ),
        "taxonomy": select_features_taxonomy_description,
    },
    name="Select downstream features from Kraken 2.",
    description=select_features_description,
)

plugin.methods.register_function(
    function=q2_annotate.kraken2.kraken2_to_mag_features,
    inputs={
        "reports": FeatureData[Kraken2Reports % Properties("mags")],
        "outputs": FeatureData[Kraken2Outputs % Properties("mags")],
    },
    parameters={
        "coverage_threshold": Float % Range(0, 100, inclusive_end=True),
        # 'lca_mode': Str % Choices(['lca', 'majority'])
    },
    outputs=[("taxonomy", FeatureData[Taxonomy])],
    input_descriptions={
        "reports": "Per-sample Kraken 2 reports.",
        "outputs": "Per-sample Kraken 2 output files.",
    },
    parameter_descriptions={
        "coverage_threshold": (
            "The minimum percent coverage required to produce a feature."
        ),
        # 'lca_mode': 'The method used to determine the LCA of a MAG using '
        #             'taxonomic assignments of its contigs. '
    },
    output_descriptions={
        "taxonomy": select_features_taxonomy_description,
    },
    name="Select downstream MAG features from Kraken 2.",
    description=select_features_description,
)

plugin.pipelines.register_function(
    function=q2_annotate.kraken2.map_taxonomy_to_contigs,
    inputs={
        "reports": SampleData[Kraken2Reports % Properties("contigs")],
        "outputs": SampleData[Kraken2Outputs % Properties("contigs")],
    },
    parameters={
        "coverage_threshold": Float % Range(0, 100, inclusive_end=True),
    },
    outputs=[
        ("feature_map", FeatureMap[TaxonomyToContigs]),
        ("taxonomy", FeatureData[Taxonomy]),
    ],
    input_descriptions={
        "reports": "Per-sample Kraken 2 reports for contigs.",
        "outputs": "Per-sample Kraken 2 output files for contigs.",
    },
    parameter_descriptions={
        "coverage_threshold": (
            "The minimum percent coverage required to produce a feature."
        ),
    },
    output_descriptions={
        "feature_map": (
            "Taxonomy assignments for contigs. Each contig ID is mapped to "
            "its full taxonomy string based on Kraken2 classifications. "
            "Unclassified contigs are assigned 'd__Unclassified'. Assignments "
            "below the provided coverage threshold will be treated as "
            "'d__Unclassified'."
        ),
        "taxonomy": select_features_taxonomy_description,
    },
    name="Map contig IDs to taxonomy strings from Kraken 2.",
    description=(
        "Maps contig IDs to their full taxonomy strings based on Kraken 2 "
        "classifications. This action processes Kraken 2 reports and outputs "
        "produced from contig sequences to create a taxonomy mapping where "
        "each contig ID is associated with its assigned taxonomy."
    ),
)

plugin.visualizers.register_function(
    function=q2_annotate.kraken2._visualize_collapsed_contigs,
    inputs={
        "table": FeatureTable[Frequency],
        "contig_map": FeatureMap[TaxonomyToContigs],
        "taxonomy": FeatureData[Taxonomy],
        "collapsed_table": FeatureTable[Frequency],
    },
    parameters={},
    input_descriptions={
        "table": (
            "Original feature table with contig IDs as feature "
            "IDs (before collapsing)."
        ),
        "contig_map": "Mapping between contig IDs and assigned taxonomy IDs.",
        "taxonomy": (
            "Optional taxonomy mapping to display taxonomy strings instead of IDs."
        ),
        "collapsed_table": "Feature table with contig IDs collapsed to taxonomy IDs.",
    },
    parameter_descriptions={},
    name="Visualize collapsed contig abundances",
    description=(
        "Generates a visualization showing histograms of contig abundance "
        "distributions per taxon per sample from the original (pre-collapsed) "
        "table."
    ),
)

plugin.pipelines.register_function(
    function=q2_annotate.kraken2.collapse_contigs,
    inputs={
        "table": FeatureTable[Frequency],
        "contig_map": FeatureMap[TaxonomyToContigs],
        "taxonomy": FeatureData[Taxonomy],
    },
    parameters={},
    outputs=[
        ("collapsed_table", FeatureTable[Frequency]),
        ("visualization", Visualization),
    ],
    input_descriptions={
        "table": (
            "Table of contig abundances per sample. Feature IDs will be "
            "replaced with taxonomy strings."
        ),
        "contig_map": ("Mapping between contig IDs and assigned taxonomy IDs."),
        "taxonomy": (
            "Optional taxonomy mapping to display taxonomy strings in the "
            "visualization instead of taxonomy IDs."
        ),
    },
    parameter_descriptions={},
    output_descriptions={
        "collapsed_table": (
            "Contig abundance table with contig IDs replaced by their "
            "taxonomy strings. Contigs with the same taxonomy assignment "
            "will have their abundances averaged."
        ),
        "visualization": (
            "Interactive visualization showing histograms of contig abundance "
            "distributions per taxon per sample."
        ),
    },
    name="Collapse the contig abundances based on taxonomy.",
    description=(
        "This action collapses contig abundances based on their taxonomy "
        "assignments. Contig abundances for contigs with the same taxonomy "
        "assignment will be averaged."
    ),
)

plugin.methods.register_function(
    function=q2_annotate.eggnog.build_custom_diamond_db,
    inputs={
        "seqs": FeatureData[ProteinSequence],
        "taxonomy": ReferenceDB[NCBITaxonomy],
    },
    input_descriptions={
        "seqs": "Protein reference database.",
        "taxonomy": ("Reference taxonomy, needed to provide taxonomy features."),
    },
    outputs=[("db", ReferenceDB[Diamond])],
    output_descriptions={"db": "DIAMOND database."},
    parameters={
        "threads": Int % Range(1, None),
        "file_buffer_size": Int % Range(1, None),
        "ignore_warnings": Bool,
        "no_parse_seqids": Bool,
    },
    parameter_descriptions={
        "threads": "Number of CPU threads.",
        "file_buffer_size": "File buffer size in bytes.",
        "ignore_warnings": "Ignore warnings.",
        "no_parse_seqids": "Print raw seqids without parsing.",
    },
    name="Create a DIAMOND formatted reference database from a FASTA input file.",
    description=(
        "Creates an artifact containing a binary DIAMOND database "
        "file (ref_db.dmnd) from a protein reference database "
        "file in FASTA format."
    ),
    citations=[citations["buchfink_sensitive_2021"]],
    examples={"Minimum working example": ex.diamond_makedb},
)

plugin.methods.register_function(
    function=q2_annotate.eggnog.fetch_eggnog_db,
    inputs={},
    parameters={},
    outputs=[("db", ReferenceDB[Eggnog])],
    output_descriptions={"db": "eggNOG annotation database."},
    name=("Fetch the databases necessary to run the eggnog-annotate action."),
    description=(
        "Downloads EggNOG reference database using the "
        "`download_eggnog_data.py` script from eggNOG. Here, this "
        "script downloads 3 files and stores them in the output "
        "artifact. At least 80 GB of storage space is required to "
        "run this action."
    ),
    citations=[citations["huerta_cepas_eggnog_2019"]],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog.fetch_diamond_db,
    inputs={},
    parameters={},
    outputs=[("db", ReferenceDB[Diamond])],
    output_descriptions={"db": "Complete Diamond reference database."},
    name=(
        "Fetch the complete Diamond database necessary to run the "
        "eggnog-diamond-search action."
    ),
    description=(
        "Downloads Diamond reference database. "
        "This action downloads 1 file (ref_db.dmnd). At least 18 GB "
        "of storage space is required to run this action."
    ),
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"],
    ],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog.fetch_eggnog_proteins,
    inputs={},
    parameters={},
    outputs=[("eggnog_proteins", ReferenceDB[EggnogProteinSequences])],
    output_descriptions={
        "eggnog_proteins": (
            "eggNOG database of protein sequences and "
            "their corresponding taxonomy information."
        )
    },
    name=("Fetch the databases necessary to run the build-eggnog-diamond-db action."),
    description=(
        "Downloads eggNOG proteome database. "
        "This script downloads 2 files "
        "(e5.proteomes.faa and e5.taxid_info.tsv) "
        "and creates and artifact with them. At least 18 GB of "
        "storage space is required to run this action."
    ),
    citations=[citations["huerta_cepas_eggnog_2019"]],
)


plugin.methods.register_function(
    function=q2_annotate.eggnog.fetch_ncbi_taxonomy,
    inputs={},
    parameters={},
    outputs=[("taxonomy", ReferenceDB[NCBITaxonomy])],
    output_descriptions={"taxonomy": "NCBI reference taxonomy."},
    name="Fetch NCBI reference taxonomy.",
    description=(
        "Downloads NCBI reference taxonomy from the NCBI FTP server. "
        "The resulting artifact is required by the "
        "build-custom-diamond-db action if one wishes to "
        "create a Diamond data base with taxonomy features. "
        "At least 30 GB of "
        "storage space is required to run this action."
    ),
    citations=[citations["NCBI"]],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog.build_eggnog_diamond_db,
    inputs={
        "eggnog_proteins": ReferenceDB[EggnogProteinSequences],
    },
    input_descriptions={
        "eggnog_proteins": (
            "eggNOG database of protein sequences and "
            "their corresponding taxonomy information "
            "(generated through the `fetch-eggnog-proteins` "
            "action)."
        )
    },
    parameters={"taxon": Int % Range(2, 1579337)},
    parameter_descriptions={"taxon": "NCBI taxon ID number."},
    outputs=[("db", ReferenceDB[Diamond])],
    output_descriptions={
        "db": "Complete Diamond reference database for the specified taxon."
    },
    name="Create a DIAMOND formatted reference database for the specified taxon.",
    description=(
        "Creates a DIAMOND database which contains the protein "
        "sequences that belong to the specified taxon."
    ),
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"],
    ],
)

plugin.pipelines.register_function(
    function=q2_annotate.eggnog.search_orthologs_diamond,
    inputs={
        "seqs": SampleData[Contigs] | SampleData[MAGs] | FeatureData[MAG],
        "db": ReferenceDB[Diamond],
    },
    parameters={"num_cpus": Int, "db_in_memory": Bool, **partition_params},
    input_descriptions={
        "seqs": "Sequences to be searched for hits using the Diamond Database",
        "db": "The filepath to an artifact containing the Diamond database",
    },
    parameter_descriptions={
        "num_cpus": ("Number of CPUs to utilize. '0' will use all available."),
        "db_in_memory": (
            "Read database into memory. The database can be very large, "
            "so this option should only be used on clusters or other "
            "machines with enough memory."
        ),
        **partition_param_descriptions,
    },
    outputs=[
        ("eggnog_hits", SampleData[Orthologs]),
        ("table", FeatureTable[Frequency]),
        ("loci", GenomeData[Loci]),
    ],
    name="Run eggNOG search using diamond aligner.",
    description=(
        "Use Diamond and eggNOG to align contig or MAG sequences "
        "against the Diamond database."
    ),
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"],
    ],
)

plugin.pipelines.register_function(
    function=q2_annotate.eggnog.search_orthologs_hmmer,
    inputs={
        "seqs": SampleData[Contigs | MAGs] | FeatureData[MAG],
        "pressed_hmm_db": ProfileHMM[PressedProtein],
        "idmap": EggnogHmmerIdmap,
        "seed_alignments": GenomeData[Proteins],
    },
    parameters={"num_cpus": Int, "db_in_memory": Bool, **partition_params},
    input_descriptions={
        "seqs": "Sequences to be searched for hits.",
        "pressed_hmm_db": "Collection of profile HMMs in binary format and indexed.",
        "idmap": "List of protein families in `pressed_hmm_db`.",
        "seed_alignments": (
            "Seed alignments for the protein families in `pressed_hmm_db`."
        ),
    },
    parameter_descriptions={
        "num_cpus": (
            "Number of CPUs to utilize per partition. '0' will use all available."
        ),
        "db_in_memory": (
            "Read database into memory. The database can be very large, "
            "so this option should only be used on clusters or other "
            "machines with enough memory."
        ),
        **partition_param_descriptions,
    },
    outputs=[
        ("eggnog_hits", SampleData[Orthologs]),
        ("table", FeatureTable[Frequency]),
        ("loci", GenomeData[Loci]),
    ],
    name="Run eggNOG search using HMMER aligner.",
    description=(
        "This method uses HMMER to find possible target sequences "
        "to annotate with eggNOG-mapper."
    ),
    citations=[
        citations["noauthor_hmmer_nodate"],
        citations["huerta_cepas_eggnog_2019"],
    ],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog._eggnog_diamond_search,
    inputs={
        "seqs": SampleData[Contigs] | SampleData[MAGs] | FeatureData[MAG],
        "db": ReferenceDB[Diamond],
    },
    parameters={"num_cpus": Int, "db_in_memory": Bool},
    input_descriptions={
        "seqs": "Sequences to be searched for ortholog hits.",
        "db": "Diamond database.",
    },
    parameter_descriptions={
        "num_cpus": ("Number of CPUs to utilize. '0' will use all available."),
        "db_in_memory": (
            "Read database into memory. The database can be very large, "
            "so this option should only be used on clusters or other "
            "machines with enough memory."
        ),
    },
    outputs=[
        ("eggnog_hits", SampleData[Orthologs]),
        ("table", FeatureTable[Frequency]),
        ("loci", GenomeData[Loci]),
    ],
    output_descriptions={
        "eggnog_hits": (
            "BLAST6-like table(s) describing the identified "
            "orthologs. One table per sample or MAG in the input."
        ),
        "table": "Feature table with counts of orthologs per sample/MAG.",
        "loci": "Loci of the identified orthologs.",
    },
    name="Run eggNOG search using Diamond aligner.",
    description=(
        "This method performs the steps by which we find our "
        "possible target sequences to annotate using the Diamond "
        "search functionality from the eggnog `emapper.py` script."
    ),
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"],
    ],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog._eggnog_hmmer_search,
    inputs={
        "seqs": SampleData[Contigs] | SampleData[MAGs] | FeatureData[MAG],
        "idmap": EggnogHmmerIdmap,
        "pressed_hmm_db": ProfileHMM[PressedProtein],
        "seed_alignments": GenomeData[Proteins],
    },
    parameters={
        "num_cpus": Int,
        "db_in_memory": Bool,
    },
    input_descriptions={
        "seqs": "Sequences to be searched for hits.",
        "idmap": "List of protein families in `pressed_hmm_db`.",
        "pressed_hmm_db": "Collection of Profile HMMs in binary format and indexed.",
        "seed_alignments": (
            "Seed alignments for the protein families in `pressed_hmm_db`."
        ),
    },
    parameter_descriptions={
        "num_cpus": (
            "Number of CPUs to utilize per partition. '0' will use all available."
        ),
        "db_in_memory": (
            "Read database into memory. The database can be very large, "
            "so this option should only be used on clusters or other "
            "machines with enough memory."
        ),
    },
    outputs=[
        ("eggnog_hits", SampleData[Orthologs]),
        ("table", FeatureTable[Frequency]),
        ("loci", GenomeData[Loci]),
    ],
    output_descriptions={
        "eggnog_hits": (
            "BLAST6-like table(s) describing the identified "
            "orthologs. One table per sample or MAG in the input."
        ),
        "table": "Feature table with counts of orthologs per sample/MAG.",
        "loci": "Loci of the identified orthologs.",
    },
    name="Run eggNOG search using HMMER aligner.",
    description=(
        "This method performs the steps by which we find our "
        "possible target sequences to annotate using the "
        "HMMER search functionality from the eggnog `emapper.py` script."
    ),
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"],
    ],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog._eggnog_feature_table,
    inputs={"seed_orthologs": SampleData[Orthologs]},
    parameters={},
    input_descriptions={
        "seed_orthologs": ("Sequence data to be turned into an eggnog feature table.")
    },
    parameter_descriptions={},
    outputs=[("table", FeatureTable[Frequency])],
    name="Create an eggnog table.",
    description="Create an eggnog table.",
)

plugin.pipelines.register_function(
    function=q2_annotate.eggnog.map_eggnog,
    inputs={
        "eggnog_hits": SampleData[Orthologs],
        "db": ReferenceDB[Eggnog],
    },
    input_descriptions={
        "eggnog_hits": ("BLAST6-like table(s) describing the identified orthologs."),
        "db": "eggNOG annotation database.",
    },
    parameters={
        "db_in_memory": Bool,
        "num_cpus": Int % Range(0, None),
        **partition_params,
    },
    parameter_descriptions={
        "db_in_memory": (
            "Read eggnog database into memory. The eggnog database is very large "
            "(>44GB), so this option should only be used on clusters or other "
            "machines with enough memory."
        ),
        "num_cpus": ("Number of CPUs to utilize. '0' will use all available."),
        **partition_param_descriptions,
    },
    outputs=[("ortholog_annotations", GenomeData[NOG])],
    output_descriptions={"ortholog_annotations": "Annotated hits."},
    name="Annotate orthologs against eggNOG database.",
    description="Apply eggnog mapper to annotate seed orthologs.",
    citations=[citations["huerta_cepas_eggnog_2019"]],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog._eggnog_annotate,
    inputs={
        "eggnog_hits": SampleData[Orthologs],
        "db": ReferenceDB[Eggnog],
    },
    parameters={"db_in_memory": Bool, "num_cpus": Int % Range(0, None)},
    parameter_descriptions={
        "db_in_memory": (
            "Read eggnog database into memory. The eggNOG database is very large "
            "(>44GB), so this option should only be used on clusters or other "
            "machines with enough memory."
        ),
        "num_cpus": ("Number of CPUs to utilize. '0' will use all available."),
    },
    outputs=[("ortholog_annotations", GenomeData[NOG])],
    name="Annotate orthologs against eggNOG database.",
    description="Apply eggnog mapper to annotate seed orthologs.",
    citations=[citations["huerta_cepas_eggnog_2019"]],
)

plugin.methods.register_function(
    function=q2_annotate.prodigal.predict_genes_prodigal,
    inputs={"seqs": FeatureData[MAG] | SampleData[MAGs] | SampleData[Contigs]},
    input_descriptions={
        "seqs": "MAGs or contigs for which one wishes to predict genes."
    },
    parameters={
        "translation_table_number": Str
        % Choices(
            [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "21",
                "22",
                "23",
                "24",
                "25",
            ]
        ),
        "mode": Str % Choices(["single", "meta"]),
        "closed": Bool,
        "no_shine_dalgarno": Bool,
        "mask": Bool,
    },
    parameter_descriptions={
        "translation_table_number": (
            "Translation table to be used to translate genes into sequences of "
            "amino acids. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/"
            "wprintgc.cgi for reference."
        ),
        "mode": (
            "Gene prediction mode. 'single' is suitable for single genome analysis "
            "(e.g., MAGs), 'meta' is suitable for metagenome analysis "
            "(e.g., contigs from mixed communities)."
        ),
        "closed": (
            "Treat sequences as complete genomes with closed ends. Use this "
            "for finished genomes where the sequences represent complete chromosomes "
            "or plasmids."
        ),
        "no_shine_dalgarno": (
            "Bypass Shine-Dalgarno trainer and use a more generic model. "
            "Useful for virus, phage, or plasmid sequences that may not follow "
            "standard prokaryotic gene patterns."
        ),
        "mask": (
            "Treat runs of N as masked sequence. Useful for assemblies that "
            "contain gap regions represented by stretches of N nucleotides."
        ),
    },
    outputs=[
        ("loci", GenomeData[Loci]),
        ("genes", GenomeData[Genes]),
        ("proteins", GenomeData[Proteins]),
    ],
    output_descriptions={
        "loci": (
            "Gene coordinates files (one per MAG or sample) listing the "
            "location of each predicted gene as well as some additional "
            "scoring information."
        ),
        "genes": (
            "Fasta files (one per MAG or sample) with the nucleotide "
            "sequences of the predicted genes."
        ),
        "proteins": (
            "Fasta files (one per MAG or sample) with the protein "
            "translation of the predicted genes."
        ),
    },
    name="Predict gene sequences from MAGs or contigs using Prodigal.",
    description=(
        "Prodigal (PROkaryotic DYnamic programming Gene-finding ALgorithm), "
        "a gene prediction algorithm designed for improved gene structure "
        "prediction, translation initiation site recognition, and reduced "
        "false positives in bacterial and archaeal genomes."
    ),
    citations=[citations["hyatt_prodigal_2010"]],
)

plugin.methods.register_function(
    function=q2_annotate.kaiju.fetch_kaiju_db,
    inputs={},
    parameters={
        "database_type": Str
        % Choices(
            [
                "nr",
                "nr_euk",
                "refseq",
                "refseq_ref",
                "refseq_nr",
                "fungi",
                "viruses",
                "plasmids",
                "progenomes",
                "rvdb",
            ]
        ),
    },
    outputs=[
        ("db", KaijuDB),
    ],
    input_descriptions={},
    parameter_descriptions={
        "database_type": (
            "Type of database to be downloaded. For more information on available "
            "types please see the list on Kaiju's web server: https://bioinformatics"
            "-centre.github.io/kaiju/downloads.html"
        ),
    },
    output_descriptions={"db": "Kaiju database."},
    name="Fetch Kaiju database.",
    description=(
        "This method fetches the latest Kaiju database from Kaiju's web server."
    ),
    citations=[citations["menzel2016"]],
)

kaiju_params = {
    "z": Int % Range(1, None),
    "a": Str % Choices(["greedy", "mem"]),
    "e": Int % Range(1, None),
    "m": Int % Range(1, None),
    "s": Int % Range(1, None),
    "evalue": Float % Range(0, 1),
    "x": Bool,
    "r": Str % Choices(["phylum", "class", "order", "family", "genus", "species"]),
    "c": Float % Range(0, 100, inclusive_start=True),
    "exp": Bool,
    "u": Bool,
}
kaiju_param_descriptions = {
    "z": "Number of threads.",
    "a": "Run mode.",
    "e": "Number of mismatches allowed in Greedy mode.",
    "m": "Minimum match length.",
    "s": "Minimum match score in Greedy mode.",
    "evalue": "Minimum E-value in Greedy mode.",
    "x": "Enable SEG low complexity filter.",
    "r": "Taxonomic rank.",
    "c": (
        "Minimum required number or fraction of reads for "
        "the taxon (except viruses) to be reported."
    ),
    "exp": (
        "Expand viruses, which are always shown as full "
        "taxon path and read counts are not summarized in "
        "higher taxonomic levels."
    ),
    "u": (
        "Do not count unclassified reads for the total reads "
        "when calculating percentages for classified reads."
    ),
}

plugin.methods.register_function(
    function=q2_annotate.kaiju._classify_kaiju,
    inputs={
        "seqs": SampleData[
            SequencesWithQuality
            | PairedEndSequencesWithQuality
            | JoinedSequencesWithQuality
        ]
        | SampleData[Contigs],
        "db": KaijuDB,
    },
    parameters=kaiju_params,
    outputs=[
        ("abundances", FeatureTable[Frequency]),
        ("taxonomy", FeatureData[Taxonomy]),
    ],
    input_descriptions={
        "seqs": "Sequences to be classified.",
        "db": "Kaiju database.",
    },
    parameter_descriptions=kaiju_param_descriptions,
    output_descriptions={
        "abundances": "Sequence abundances.",
        "taxonomy": "Linked taxonomy.",
    },
    name="Classify sequences using Kaiju.",
    description=(
        "This method uses Kaiju to perform taxonomic classification "
        "of DNA sequence reads or contigs."
    ),
    citations=[citations["menzel2016"]],
)

plugin.pipelines.register_function(
    function=q2_annotate.kaiju.classify_kaiju,
    inputs={
        "seqs": SampleData[
            SequencesWithQuality
            | PairedEndSequencesWithQuality
            | JoinedSequencesWithQuality
        ]
        | SampleData[Contigs],
        "db": KaijuDB,
    },
    parameters={**kaiju_params, **partition_params},
    outputs=[
        ("abundances", FeatureTable[Frequency]),
        ("taxonomy", FeatureData[Taxonomy]),
    ],
    input_descriptions={
        "seqs": "Sequences to be classified.",
        "db": "Kaiju database.",
    },
    parameter_descriptions={**kaiju_param_descriptions, **partition_param_descriptions},
    output_descriptions={
        "abundances": "Sequence abundances.",
        "taxonomy": "Linked taxonomy.",
    },
    name="Classify sequences using Kaiju.",
    description="This method uses Kaiju to perform taxonomic classification.",
    citations=[citations["menzel2016"]],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog.fetch_eggnog_hmmer_db,
    inputs={},
    parameters={"taxon_id": Int % Range(2, None)},
    parameter_descriptions={"taxon_id": "Taxon ID number."},
    outputs=[
        ("idmap", EggnogHmmerIdmap % Properties("eggnog")),
        ("hmm_db", ProfileHMM[MultipleProtein] % Properties("eggnog")),
        ("pressed_hmm_db", ProfileHMM[PressedProtein] % Properties("eggnog")),
        ("seed_alignments", GenomeData[Proteins] % Properties("eggnog")),
    ],
    output_descriptions={
        "idmap": "List of protein families in `hmm_db`.",
        "hmm_db": "Collection of Profile HMMs.",
        "pressed_hmm_db": ("Collection of Profile HMMs in binary format and indexed."),
        "seed_alignments": ("Seed alignments for the protein families in `hmm_db`."),
    },
    name=(
        "Fetch the taxon specific database necessary to run the "
        "eggnog-hmmer-search action."
    ),
    description="Downloads Profile HMM database for the specified taxon.",
    citations=[
        citations["huerta_cepas_eggnog_2019"],
        citations["noauthor_hmmer_nodate"],
    ],
)

I_reads, O_reads = TypeMap(
    {
        SampleData[SequencesWithQuality]: SampleData[SequencesWithQuality],
        SampleData[PairedEndSequencesWithQuality]: SampleData[
            PairedEndSequencesWithQuality
        ],
    }
)

plugin.pipelines.register_function(
    function=q2_annotate.filtering.construct_human_pangenome_index,
    inputs={},
    parameters={"threads": Threads},
    outputs=[("index", Bowtie2Index)],
    input_descriptions={},
    parameter_descriptions={
        "threads": "Number of threads to use when building the index."
    },
    output_descriptions={"index": "Generated combined human reference index."},
    name="Construct the human pangenome index.",
    description=(
        "This method generates a Bowtie2 index for the combined human "
        "GRCh38 reference genome and the draft human pangenome."
    ),
    citations=[],
)

plugin.pipelines.register_function(
    function=q2_annotate.filtering.filter_reads_human_pangenome,
    inputs={"reads": I_reads, "index": Bowtie2Index},
    parameters={
        "threads": Threads,
        **{
            k: v
            for (k, v) in filter_parameters.items()
            if k not in ["exclude_seqs", "n_threads"]
        },
    },
    outputs=[("filtered_reads", O_reads), ("reference_index", Bowtie2Index)],
    input_descriptions={
        "reads": "Reads to be filtered against the human genome.",
        "index": (
            "Bowtie2 index of the reference human genome. If not "
            "provided, an index combined from the reference GRCh38 "
            "human genome and the human pangenome will be generated."
        ),
    },
    parameter_descriptions={
        "threads": "Number of threads to use for indexing and read filtering.",
        **{
            k: v
            for (k, v) in filter_parameter_descriptions.items()
            if k not in ["exclude_seqs", "n_threads"]
        },
    },
    output_descriptions={
        "filtered_reads": ("Original reads without the contaminating human reads."),
        "reference_index": (
            "Generated combined human reference index. If an "
            "index was provided as an input, it will be "
            "returned here instead."
        ),
    },
    name="Remove contaminating human reads.",
    description=(
        "Generates a Bowtie2 index fo the combined human "
        "GRCh38 reference genome and the draft human pangenome, and"
        "uses that index to remove the contaminating human reads from "
        "the reads provided as input."
    ),
    citations=[],
)

plugin.methods.register_function(
    function=q2_annotate.eggnog.annotation.extract_annotations,
    inputs={
        "ortholog_annotations": GenomeData[NOG],
    },
    parameters={
        "annotation": Str
        % Choices(
            [
                "cog",
                "caz",
                "kegg_ko",
                "kegg_pathway",
                "kegg_reaction",
                "kegg_module",
                "brite",
                "ec",
            ]
        ),
        "max_evalue": Float % Range(0, None),
        "min_score": Float % Range(0, None),
    },
    outputs=[("annotation_frequency", FeatureTable[Frequency])],
    input_descriptions={"ortholog_annotations": "Ortholog annotations."},
    parameter_descriptions={"annotation": "Annotation to extract."},
    output_descriptions={
        "annotation_frequency": ("Feature table with frequency of each annotation."),
    },
    name="Extract annotation frequencies from all annotations.",
    description=(
        "This method extract a specific annotation from the table "
        "generated by EggNOG and calculates its frequencies across "
        "all MAGs."
    ),
    citations=[],
)

TMR = TypeMatch(
    [
        SampleData[Kraken2Reports % Properties("reads")],
        SampleData[Kraken2Reports % Properties("contigs")],
        SampleData[Kraken2Reports % Properties("mags")],
        FeatureData[Kraken2Reports % Properties("mags")],
    ]
)
TMO = TypeMatch(
    [
        SampleData[Kraken2Outputs % Properties("reads")],
        SampleData[Kraken2Outputs % Properties("contigs")],
        SampleData[Kraken2Outputs % Properties("mags")],
        FeatureData[Kraken2Outputs % Properties("mags")],
    ]
)


plugin.methods.register_function(
    function=q2_annotate.kraken2._filter_kraken2_reports_by_abundance,
    inputs={
        "reports": TMR,
    },
    parameters={
        "abundance_threshold": Float % Range(0, 1, inclusive_end=True),
        "remove_empty": Bool,
    },
    outputs=[("filtered_reports", TMR)],
    input_descriptions={
        "reports": "The kraken2 reports to filter by relative abundance."
    },
    parameter_descriptions={
        "abundance_threshold": (
            "A proportion between 0 and 1 representing the minimum relative "
            "abundance (by *classified* read count) that a taxon must have to "
            "be retained in the filtered report."
        ),
        "remove_empty": (
            "If True, reports with only unclassified reads "
            "remaining will be removed from the filtered data."
        ),
    },
    output_descriptions={
        "filtered_reports": "The relative abundance-filtered kraken2 reports"
    },
    name="Filter kraken2 reports by relative abundance.",
    description=(
        "Filters kraken2 reports on a per-taxon basis by relative abundance "
        "(relative frequency). Useful for removing suspected spurious "
        "classifications."
    ),
    citations=[],
)

filter_kraken2_results_param_desc = {
    "metadata": (
        "Metadata indicating which IDs to filter. The optional "
        "`where` parameter may be used to filter IDs based on "
        "specified conditions in the metadata. The optional "
        "`exclude_ids` parameter may be used to exclude the IDs "
        "specified in the metadata from the filter."
    ),
    "where": (
        "Optional SQLite WHERE clause specifying metadata criteria that "
        "must be met to be included in the filtered data. If not "
        "provided, all IDs in `metadata` that are also in the data will "
        "be retained."
    ),
    "exclude_ids": (
        "If True, the samples selected by the `metadata` and "
        "optional `where` parameter will be excluded from the "
        "filtered data."
    ),
    "remove_empty": (
        "If True, reports with only unclassified reads will be "
        "removed from the filtered data. Reports containing "
        "sequences classified only as root will also be removed."
    ),
}

plugin.methods.register_function(
    function=q2_annotate.kraken2._filter_kraken2_results_by_metadata,
    inputs={
        "reports": TMR,
        "outputs": TMO,
    },
    parameters={
        "metadata": Metadata,
        "where": Str,
        "exclude_ids": Bool,
        "remove_empty": Bool,
    },
    outputs={
        "filtered_reports": TMR,
        "filtered_outputs": TMO,
    },
    input_descriptions={
        "reports": "The Kraken reports to filter.",
        "outputs": "The Kraken outputs to filter.",
    },
    parameter_descriptions=filter_kraken2_results_param_desc,
    name="Filter Kraken2 reports and outputs.",
    description=(
        "Filter Kraken2 reports and outputs based on metadata or remove "
        "reports with 100% unclassified reads."
    ),
)

KRAKEN2_REPORTS_T = TypeMatch(
    [
        SampleData[Kraken2Reports % Properties("reads")],
        SampleData[Kraken2Reports % Properties("contigs")],
        SampleData[Kraken2Reports % Properties("mags")],
        FeatureData[Kraken2Reports % Properties("mags")],
    ]
)
KRAKEN2_OUTPUTS_T = TypeMatch(
    [
        SampleData[Kraken2Outputs % Properties("reads")],
        SampleData[Kraken2Outputs % Properties("contigs")],
        SampleData[Kraken2Outputs % Properties("mags")],
        FeatureData[Kraken2Outputs % Properties("mags")],
    ]
)
plugin.methods.register_function(
    function=q2_annotate.kraken2._merge_kraken2_results,
    inputs={"reports": List[KRAKEN2_REPORTS_T], "outputs": List[KRAKEN2_OUTPUTS_T]},
    parameters={},
    outputs={"merged_reports": KRAKEN2_REPORTS_T, "merged_outputs": KRAKEN2_OUTPUTS_T},
    input_descriptions={
        "reports": (
            "The kraken2 reports to merge. Only reports with the same sample "
            "ID are merged into one report."
        ),
        "outputs": (
            "The kraken2 outputs to merge. Only outputs with the same sample "
            "ID are merged into one output."
        ),
    },
    parameter_descriptions={},
    output_descriptions={
        "merged_reports": "The merged kraken2 reports.",
        "merged_outputs": "The merged kraken2 outputs.",
    },
    name="Merge kraken2 reports and outputs.",
    description=(
        "Merge multiple kraken2 reports and outputs such that the results "
        "contain a union of the samples represented in the inputs. "
        "If sample IDs overlap across the inputs, these reports and outputs "
        "will be processed into a single report or output per sample ID."
    ),
)

plugin.methods.register_function(
    function=q2_annotate.kraken2._align_outputs_with_reports,
    inputs={
        "reports": TMR,
        "outputs": TMO,
    },
    parameters={},
    outputs=[("aligned_outputs", TMO)],
    input_descriptions={
        "reports": "The filtered kraken2 reports.",
        "outputs": "The kraken2 outputs to align with the filtered reports.",
    },
    output_descriptions={
        "aligned_outputs": "The report-aligned filtered kraken2 outputs."
    },
    name="Align unfiltered kraken2 outputs with filtered kraken2 reports.",
    description="",
)

plugin.pipelines.register_function(
    function=q2_annotate.kraken2.filter_kraken2_results,
    inputs={
        "reports": TMR,
        "outputs": TMO,
    },
    parameters={
        "metadata": Metadata,
        "where": Str,
        "exclude_ids": Bool,
        "remove_empty": Bool,
        "abundance_threshold": Float % Range(0, 1, inclusive_end=True),
    },
    outputs={"filtered_reports": TMR, "filtered_outputs": TMO},
    input_descriptions={
        "reports": "The kraken2 reports to filter.",
        "outputs": "The kraken2 outputs to filter.",
    },
    parameter_descriptions={
        "metadata": filter_kraken2_results_param_desc["metadata"],
        "where": filter_kraken2_results_param_desc["where"],
        "exclude_ids": filter_kraken2_results_param_desc["exclude_ids"],
        "remove_empty": filter_kraken2_results_param_desc["remove_empty"],
        "abundance_threshold": (
            "A proportion between 0 and 1 representing the minimum relative "
            "abundance (by *classified* read count) that a taxon must have to "
            "be retained in the filtered report. If a taxon is filtered from "
            "the report, its associated read counts are removed entirely from "
            "the report (i.e., the subtraction of those counts is propagated "
            "to parent taxonomic groupings)."
        ),
    },
    output_descriptions={
        "filtered_reports": "The filtered kraken2 reports.",
        "filtered_outputs": "The filtered kraken2 outputs.",
    },
    name="Filter kraken2 reports and outputs by metadata and abundance.",
    description=(
        "Filter kraken2 reports and outputs by sample metadata, and/or filter "
        "classified taxa by relative abundance."
    ),
)

filter_reads_kraken2_params = {
    "taxonomy": Str,
    "include_descendants": Bool,
    "contains": Bool,
    "exclude": Bool,
}
filter_reads_kraken2_param_desc = {
    "taxonomy": (
        "Taxonomy query used for read filtering. Can be a Kraken2 taxon name "
        '(for example, "Bacteria") or a taxon ID (for example, "2").'
    ),
    "include_descendants": (
        "If True, include all descendant taxa of each matching taxon."
    ),
    "contains": (
        "If True, match taxon names using case-insensitive substring matching. "
        "If False, use exact case-insensitive matching."
    ),
    "exclude": (
        "If False, retain reads that match the taxonomy query. "
        "If True, discard matching reads and retain the rest."
    ),
}
filter_reads_kraken2_input_desc = {
    "reads": (
        "The original reads that were classified by Kraken2. "
        "The sample IDs and read headers must match those used to "
        "generate `reports` and `outputs`."
    ),
    "reports": (
        "Kraken2 reports generated from `reads`. Used to identify matching "
        "taxa and (optionally) all descendant taxa."
    ),
    "outputs": (
        "Kraken2 per-read outputs generated from `reads`. Used to map matched "
        "taxa to read IDs that will be filtered."
    ),
}

plugin.methods.register_function(
    function=q2_annotate.kraken2._filter_reads_kraken2,
    inputs={
        "reads": I_reads,
        "reports": SampleData[Kraken2Reports % Properties("reads")],
        "outputs": SampleData[Kraken2Outputs % Properties("reads")],
    },
    parameters=filter_reads_kraken2_params,
    outputs=[("filtered_reads", O_reads)],
    input_descriptions=filter_reads_kraken2_input_desc,
    parameter_descriptions=filter_reads_kraken2_param_desc,
    output_descriptions={
        "filtered_reads": "Reads filtered according to Kraken2 taxonomy matches."
    },
    name="Filter Kraken2-classified reads by taxonomy.",
    description=(
        "Filter single-end or paired-end reads by matching Kraken2-assigned "
        "taxonomy, with optional descendant expansion and inverse filtering."
    ),
)

plugin.pipelines.register_function(
    function=q2_annotate.kraken2.filter_reads_kraken2,
    inputs={
        "reads": I_reads,
        "reports": SampleData[Kraken2Reports % Properties("reads")],
        "outputs": SampleData[Kraken2Outputs % Properties("reads")],
    },
    parameters={**filter_reads_kraken2_params, **partition_params},
    outputs=[("filtered_reads", O_reads)],
    input_descriptions=filter_reads_kraken2_input_desc,
    parameter_descriptions={
        **filter_reads_kraken2_param_desc,
        **partition_param_descriptions,
    },
    output_descriptions={
        "filtered_reads": "Reads filtered according to Kraken2 taxonomy matches."
    },
    name="Filter Kraken2-classified reads by taxonomy.",
    description=(
        "Filter single-end or paired-end reads by matching Kraken2-assigned "
        "taxonomy, with optional descendant expansion and inverse filtering."
    ),
)


plugin.register_formats(EggnogHmmerIdmapFileFmt, EggnogHmmerIdmapDirectoryFmt)
plugin.register_semantic_types(EggnogHmmerIdmap)
plugin.register_semantic_type_to_format(EggnogHmmerIdmap, EggnogHmmerIdmapDirectoryFmt)
