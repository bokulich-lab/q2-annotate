# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import altair as alt
import pandas as pd

alt.data_transformers.disable_max_rows()


def _draw_detailed_plots(
    df: pd.DataFrame,
    height: int = None,
    label_font_size: int = None,
    title_font_size: int = None,
    assembly_metric: str = "scaffold_n50",
    show_mag_labels: bool = False,
) -> dict:
    """
    Draws a horizontal normalized bar plot for every sample for which BUSCO was
    run. Each barplot shows the BUSCO results for each of the MAGs in the
    sample. The plots for all samples are drawn in one composite plot which
    is then returned as a dictionary for rendering (but casted to a string).

    Args:
        df (pd.DataFrame): tabular batch summary for all samples
        height (int): height of each bar in the plot
        label_font_size (int): size of the labels in plot
        title_font_size (int): size of titles in plot
    Output:
        Output plot in dictionary from casted to a string.
    """
    # Prepare data
    busco_plot_data = pd.melt(
        df,
        id_vars=["sample_id", "mag_id", "dataset", "n_markers"],
        value_vars=["single", "duplicated", "fragmented", "missing"],
        value_name="BUSCO_percentage",
        var_name="category",
    )

    secondary_plot_data = df[
        [
            "sample_id",
            "mag_id",
            "scaffold_n50",
            "contigs_n50",
            "percent_gaps",
            "scaffolds",
        ]
    ]

    # Specify order
    mapping = {"single": 1, "duplicated": 2, "fragmented": 3, "missing": 4}
    busco_plot_data["order"] = busco_plot_data["category"].map(mapping)

    # Estimate fraction of sequences in each BUSCO category
    busco_plot_data["frac_markers"] = (
        "~"
        + round(
            busco_plot_data["BUSCO_percentage"] * busco_plot_data["n_markers"] / 100
        )
        .map(int)
        .map(str)
        + "/"
        + busco_plot_data["n_markers"].map(str)
    )

    # Individual sample approach - no faceting, just simple plots for one sample
    # Create the BUSCO plot (no faceting since we're handling one sample at a time)
    # For hconcat to work properly with container width, we need to avoid "container"
    # on individual charts and use explicit widths that will scale proportionally
    busco_plot = (
        alt.Chart(busco_plot_data)
        .mark_bar()
        .encode(
            x=alt.X("sum(BUSCO_percentage)", stack="normalize", title="BUSCO fraction"),
            y=alt.Y("mag_id", axis=alt.Axis(titleFontSize=0, labels=show_mag_labels)),
            color=alt.Color(
                "category",
                scale=alt.Scale(
                    domain=["single", "duplicated", "fragmented", "missing"],
                    range=[
                        "#4A90A4",
                        "#F5A623",
                        "#D2691E",
                        "#A0522D",
                    ],  # Soft Teal, Soft Gold, Soft Orange, Soft Brown
                ),
                legend=None,  # We'll create HTML legend
            ),
            order=alt.Order("order", sort="ascending"),
            tooltip=[
                alt.Tooltip("sample_id", title="Sample ID"),
                alt.Tooltip("mag_id", title="MAG ID"),
                alt.Tooltip("dataset", title="Lineage dataset"),
                alt.Tooltip(
                    "frac_markers", title="Approx. number of markers in this category"
                ),
                alt.Tooltip("BUSCO_percentage", title="% BUSCOs"),
            ],
            opacity=alt.value(0.85),
        )
        # Don't set width here - let hconcat handle distribution
        .properties(height={"step": height})
    )

    # Create the assembly statistics plot (no faceting) - using dynamic metric
    metric_titles = {
        "scaffold_n50": "Scaffold N50 (bp)",
        "contigs_n50": "Contig N50 (bp)",
        "percent_gaps": "Percent Gaps (%)",
        "scaffolds": "Number of Scaffolds",
    }

    secondary_plot = (
        alt.Chart(secondary_plot_data)
        .mark_bar()
        .encode(
            x=alt.X(f"{assembly_metric}:Q").title(
                metric_titles.get(assembly_metric, "Assembly Metric")
            ),
            y=alt.Y("mag_id", axis=None),
            tooltip=[
                alt.Tooltip(
                    f"{assembly_metric}:Q",
                    title=metric_titles.get(assembly_metric, "Value"),
                )
            ],
            opacity=alt.value(0.85),
        )
        # Don't set width here - let hconcat handle distribution
        .properties(height={"step": height})
    )

    # Concatenate plots horizontally
    # Use spacing to control gap between plots
    output_plot = (
        alt.hconcat(busco_plot, secondary_plot, spacing=3)
        .configure_axis(labelFontSize=label_font_size, titleFontSize=title_font_size)
        .configure_legend(labelFontSize=label_font_size, titleFontSize=title_font_size)
        .configure_header(labelFontSize=label_font_size, titleFontSize=title_font_size)
    )

    spec = output_plot.to_dict()

    # Use responsive width with autosize - Vega will handle resizing
    # Set a base width that will scale proportionally
    spec["width"] = "container"  # Base width for both plots (600 + 200)
    # spec["autosize"] = {"type": "fit", "contains": "padding", "resize": True}

    # Set explicit widths on hconcat children to maintain 3:1 ratio
    # These will scale proportionally when autosize resizes
    # if "hconcat" in spec and len(spec["hconcat"]) == 2:
    #     # 600:200 ratio maintains 3:1 proportion (total 800)
    #     spec["hconcat"][0]["width"] = 600
    #     spec["hconcat"][1]["width"] = 200

    return spec
