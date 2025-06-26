# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import List

import altair as alt
import pandas as pd

alt.data_transformers.disable_max_rows()

LABEL_FONT_SIZE = 12
TITLE_FONT_SIZE = 15
PLOT_DIM = 260


def _draw_horizontal_histograms(data: pd.DataFrame, columns: List[str]):
    data = pd.melt(
        data,
        id_vars=["sample_id", "mag_id", "dataset", "n_markers"],
        value_vars=columns,
        value_name="metric",
        var_name="category",
    )

    charts = []
    for i, category in enumerate(columns):
        x_title = category.replace("_", " ").capitalize()
        y_title = "MAG count" if i == 0 else None
        chart = (
            alt.Chart(data[data["category"] == category])
            .mark_bar()
            .encode(
                x=alt.X("metric:Q", bin=True, title=x_title),
                y=alt.Y("count()", title=y_title),
                color=alt.value("steelblue"),
            )
            .properties(width=PLOT_DIM, height=PLOT_DIM)
        )
        charts.append(chart)

    chart = alt.hconcat(*charts)

    return chart


def _draw_marker_summary_histograms(data: pd.DataFrame) -> dict:
    """
    Draws summary histograms for the BUSCO marker results of all samples.

    Returns:
        dict: Dictionary containing the Vega spec.
    """
    cols = [
        ["single", "duplicated", "fragmented", "missing", "completeness"],
        ["contamination", "scaffolds", "contigs_n50", "scaffold_n50", "length"],
    ]

    if not ("completeness" in data.columns and "contamination" in data.columns):
        cols[0].remove("completeness")
        cols[1].remove("contamination")

    chart = _draw_horizontal_histograms(data, columns=cols[0])
    chart2 = _draw_horizontal_histograms(data, columns=cols[1])

    chart = (
        alt.vconcat(chart, chart2)
        .configure_axis(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_legend(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_header(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
    )

    return chart.to_dict()


def _draw_completeness_vs_contamination(data: pd.DataFrame):
    """
    Draws scatterplot of completeness vs. contamination. The user can choose to
    display all MAGs or choose one with a dropdown menu.

    Returns:
        dict: Dictionary containing the Vega spec.
    """
    color_field = "sample_id" if data["sample_id"].notnull().all() else "mag_id"
    color_title = "Sample ID" if color_field == "sample_id" else "MAG ID"

    tooltip = [
        f"{col}:Q" if pd.api.types.is_numeric_dtype(data[col]) else f"{col}:N"
        for col in data.columns
    ]

    chart = alt.Chart(data)

    unique_ids = sorted(data[color_field].dropna().unique().tolist())
    selection = alt.param(
        name="selected_id",
        bind=alt.binding_select(options=["All"] + unique_ids, name=f"{color_title}: "),
        value="All",
    )

    chart = chart.transform_filter(
        f"(selected_id == 'All') || (datum.{color_field} == selected_id)"
    ).add_params(selection)

    chart = (
        chart.mark_circle(size=60)
        .encode(
            x=alt.X(
                "completeness:Q", title="Completeness", scale=alt.Scale(domain=[0, 100])
            ),
            y=alt.Y(
                "contamination:Q",
                title="Contamination",
                scale=alt.Scale(domain=[0, 100]),
            ),
            color=alt.Color(
                f"{color_field}:N", title=color_title, scale=alt.Scale(scheme="viridis")
            ),
            tooltip=tooltip,
        )
        .properties(width=600, height=600)
        .configure_axis(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_legend(
            labelFontSize=LABEL_FONT_SIZE,
            titleFontSize=TITLE_FONT_SIZE,
            labelLimit=1000,
        )
        .configure_header(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .interactive()
    )

    return chart.to_dict()


def _draw_selectable_summary_histograms(data: pd.DataFrame) -> dict:
    """
    Draws summary histograms for the MAG assembly metrics where users
    can indicate which metric and for which sample they want to see.

    Returns:
        dict: Dictionary containing the Vega spec.
    """
    metrics = [
        "single",
        "duplicated",
        "fragmented",
        "missing",
        "completeness",
        "contamination",
        "scaffolds",
        "contigs_n50",
        "scaffold_n50",
        "length",
    ]

    if not ("completeness" in data.columns and "contamination" in data.columns):
        metrics.remove("completeness")
        metrics.remove("contamination")

    data = pd.melt(
        data,
        id_vars=["sample_id", "mag_id", "dataset", "n_markers"],
        value_vars=metrics,
        value_name="metric",
        var_name="category",
    )

    # Create the dropdown selection with all possible metrics
    selection_metrics = alt.selection_point(
        fields=["category"],
        bind=alt.binding_select(options=metrics),
        name="select_metric",
        value="single",
    )

    # Create the sample search box
    samples = data["sample_id"].unique().tolist()
    selection_box = alt.param(
        value=samples[0],
        name="select_sample",
        bind=alt.binding(
            input="search",
            placeholder="Sample IDs",
            name="select_sample",
        ),
    )

    # Create the chart
    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x=alt.X("metric:Q", bin=True, title=None),
            y=alt.Y("count()", title="MAG count"),
            color=alt.value("steelblue"),
        )
        .add_params(selection_metrics, selection_box)
        .transform_filter(
            "datum.sample_id == trim(select_sample) "
            "& datum.category == select_metric.category"
        )
        .configure_axis(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_legend(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_header(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .properties(width=PLOT_DIM, height=PLOT_DIM)
    )

    return chart.to_dict()
