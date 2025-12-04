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
PLOT_DIM = 225
BASE_COLOR = "#2d718e"
VIOLIN_COLOR = "#29af7f"


def _draw_horizontal_histogram_spec(data: pd.DataFrame, category: str) -> dict:
    """
    Build a single, responsive Vega-Lite spec for one metric/category histogram.
    Includes a shared 'selected_id' param to allow external filtering by sample ID.
    """
    x_title = category.replace("_", " ").capitalize()

    selection = alt.param(name="selected_id", value="All")

    chart = (
        alt.Chart(data[data["category"] == category])
        .mark_bar()
        .encode(
            x=alt.X("metric:Q", bin=True, title=x_title),
            y=alt.Y("count()", title="MAG count"),
            color=alt.value(BASE_COLOR),
        )
        .add_params(selection)
        .transform_filter("(selected_id == 'All') || (datum.sample_id == selected_id)")
        .properties(height=PLOT_DIM, width="container")
        .configure_axis(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_header(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
    )

    return chart.to_dict()


def _draw_marker_summary_histograms(data: pd.DataFrame) -> list:
    """
    Build responsive, per-metric histogram specs for BUSCO markers and
    assembly metrics. Returns a list of Vega-Lite specs (dicts).
    """
    cols = [
        ["single", "duplicated", "fragmented", "missing", "completeness"],
        ["contamination", "contigs_n50", "length"],
    ]

    if not ("completeness" in data.columns and "contamination" in data.columns):
        cols[0].remove("completeness")
        cols[1].remove("contamination")

    melted = pd.melt(
        data,
        id_vars=["sample_id", "mag_id", "dataset", "n_markers"],
        value_vars=[*cols[0], *cols[1]],
        value_name="metric",
        var_name="category",
    )

    specs: List[dict] = []
    for category in [*cols[0], *cols[1]]:
        specs.append(_draw_horizontal_histogram_spec(melted, category))

    return specs


def _draw_completeness_vs_contamination(data: pd.DataFrame):
    """
    Draws scatterplot of completeness vs. contamination. Filtering is controlled by
    an external 'selected_id' param (no bound UI here) so that one dropdown can drive
    multiple charts on the page.

    Returns:
        dict: Dictionary containing the Vega spec.
    """
    color_field = "sample_id" if data["sample_id"].notnull().all() else "mag_id"
    color_title = "Sample ID" if color_field == "sample_id" else "MAG ID"

    tooltip = [
        f"{col}:Q" if pd.api.types.is_numeric_dtype(data[col]) else f"{col}:N"
        for col in data.columns
    ]

    # Calculate data-driven axis upper bounds (10% above max, capped at 110)
    max_comp = (
        0
        if "completeness" not in data.columns
        else pd.to_numeric(data["completeness"], errors="coerce").max(skipna=True)
    )
    max_cont = (
        0
        if "contamination" not in data.columns
        else pd.to_numeric(data["contamination"], errors="coerce").max(skipna=True)
    )
    max_comp = 0 if pd.isna(max_comp) else float(max_comp)
    max_cont = 0 if pd.isna(max_cont) else float(max_cont)
    upper_x = min(110.0, round(max_comp * 1.1, 1))
    upper_y = min(110.0, round(max_cont * 1.1, 1))
    # Ensure at least [0,5] if all zeros
    upper_x = max(5.0, upper_x)
    upper_y = max(5.0, upper_y)

    chart = alt.Chart(data)

    selection = alt.param(name="selected_id", value="All")

    chart = chart.transform_filter(
        f"(selected_id == 'All') || (datum.{color_field} == selected_id)"
    ).add_params(selection)

    chart = (
        chart.mark_circle(size=100)
        .encode(
            x=alt.X(
                "completeness:Q",
                title="Completeness",
                scale=alt.Scale(domain=[0, upper_x]),
            ),
            y=alt.Y(
                "contamination:Q",
                title="Contamination",
                scale=alt.Scale(domain=[0, upper_y]),
            ),
            color=alt.Color(
                f"{color_field}:N",
                title=color_title,
                scale=alt.Scale(scheme="viridis"),
                legend=alt.Legend(orient="right"),
            ),
            tooltip=tooltip,
        )
        .configure_axis(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_legend(
            labelFontSize=LABEL_FONT_SIZE,
            titleFontSize=TITLE_FONT_SIZE,
            labelLimit=1000,
        )
        .configure_header(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .properties(height=400, width="container")
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
        "contigs_n50",
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
            color=alt.value(BASE_COLOR),
        )
        .add_params(selection_metrics, selection_box)
        .transform_filter(
            "datum.sample_id == trim(select_sample) "
            "& datum.category == select_metric.category"
        )
        .configure_axis(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_header(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .properties(height=300, width="container")
    )

    return chart.to_dict()


def _draw_box_whiskers_plots(data: pd.DataFrame) -> list:
    """
    Draws individual box-whiskers plot specs for BUSCO marker categories and
    assembly metrics. Each category gets its own responsive chart spec with
    individual y-axis scale and data points overlaid. Sample selection highlights
    the selected sample across all box plots.

    Returns:
        list: List of Vega-Lite specs (dicts) for individual box plots.
    """
    metrics = [
        "single",
        "duplicated",
        "fragmented",
        "missing",
        "completeness",
        "contamination",
        "contigs_n50",
        "length",
    ]

    # Remove metrics that don't exist in the data
    if not ("completeness" in data.columns and "contamination" in data.columns):
        metrics.remove("completeness")
        metrics.remove("contamination")

    # Create selection parameter for sample filtering
    selection = alt.param(name="selected_id", value="All")

    specs: List[dict] = []
    for metric in metrics:
        # Determine the appropriate y-axis format based on metric type
        if metric in ["single", "duplicated", "fragmented", "missing"]:
            # BUSCO marker metrics: no decimals, no scientific notation
            y_format = ".0f"
        else:
            # Assembly metrics: use abbreviated format (K, M, etc.)
            y_format = ".2s"

        # Create base chart with box plot for this specific metric
        # (always show all data)
        base_chart = (
            alt.Chart(data)
            .mark_boxplot(size=40)
            .encode(
                y=alt.Y(
                    f"{metric}:Q",
                    title=metric,
                    scale=alt.Scale(zero=False),
                    axis=alt.Axis(format=y_format),
                ),
                color=alt.value(BASE_COLOR),
                opacity=alt.value(0.3),
                stroke=alt.value("black"),
                strokeWidth=alt.value(0.75),
            )
            .add_params(selection)
        )

        # Create overlay chart with individual points for this specific metric
        # (show all points, highlight selected)
        points_chart = (
            alt.Chart(data)
            .mark_circle(size=50, opacity=0.4, stroke=None)
            .encode(
                y=alt.Y(
                    f"{metric}:Q",
                    title=metric,
                    scale=alt.Scale(zero=False),
                    axis=alt.Axis(format=y_format),
                ),
                color=alt.condition(
                    alt.datum.sample_id == alt.param("selected_id"),
                    alt.value("red"),
                    alt.value("gray"),
                ),
                opacity=alt.condition(
                    alt.datum.sample_id == alt.param("selected_id"),
                    alt.value(0.9),
                    alt.value(0.4),
                ),
                tooltip=[
                    alt.Tooltip("sample_id:N", title="Sample ID"),
                    alt.Tooltip("mag_id:N", title="MAG ID"),
                    alt.Tooltip(f"{metric}:Q", title="Value", format=".2f"),
                ],
            )
            .add_params(selection)
            # No filter - always show all points
        )

        # Combine box plot and points
        chart = (base_chart + points_chart).resolve_scale(y="shared")

        # Configure the chart
        chart = (
            chart.properties(height=300, width="container")
            .configure_axis(
                labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE, grid=False
            )
            .configure_header(
                labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE
            )
        )
        specs.append(chart.to_dict())

    return specs


def _draw_selectable_unbinned_histograms(data: pd.DataFrame) -> dict:
    """
    Draws summary histograms for unbinned contigs (percentage only) across samples.
    Returns a single responsive Vega spec with counts on the y-axis.
    """
    # Keep only one row per sample
    data = data.drop_duplicates(subset=["sample_id"])

    chart = (
        alt.Chart(data)
        .mark_bar()
        .encode(
            x=alt.X(
                "unbinned_contigs_count:Q", bin=True, title="Unbinned contig count"
            ),
            y=alt.Y(
                "count()",
                title="Sample count",
                axis=alt.Axis(format=".0f", tickMinStep=1),
                scale=alt.Scale(zero=True),
            ),
            color=alt.value(BASE_COLOR),
        )
        .configure_axis(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .configure_header(labelFontSize=LABEL_FONT_SIZE, titleFontSize=TITLE_FONT_SIZE)
        .properties(height=400, width="container")
    )

    return chart.to_dict()
