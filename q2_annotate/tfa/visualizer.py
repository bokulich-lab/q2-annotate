# -----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
"""Interactive exploration of sample-resolved TFA loads."""

import json
import shutil
from collections import Counter
from importlib import resources
from pathlib import Path

import biom
import numpy as np
import pandas as pd
from qiime2 import Metadata

from ._mapping import _align_gene_taxonomy


def _taxon_labels(pairs: list[tuple[str, str]], taxonomy: pd.DataFrame):
    ids = list(dict.fromkeys(taxon_id for taxon_id, _ in pairs))
    full = {}
    for taxon_id in ids:
        full[taxon_id] = str(taxonomy.at[taxon_id, "Taxon"]).strip()

    duplicates = Counter(full.values())
    full = {
        taxon_id: f"{label} [{taxon_id}]" if duplicates[label] > 1 else label
        for taxon_id, label in full.items()
    }
    short = {
        taxon_id: label.rsplit(";", 1)[-1].strip() for taxon_id, label in full.items()
    }
    short_duplicates = Counter(short.values())
    short = {
        taxon_id: f"{label} [{taxon_id}]" if short_duplicates[label] > 1 else label
        for taxon_id, label in short.items()
    }
    return full, short


def _sample_summaries(values: np.ndarray, sample_count: int) -> dict:
    total = float(values.sum())
    if not sample_count:
        return {"total": total, "mean": 0.0, "median": 0.0}

    # Nonnegative loads place implicit sparse zeros before stored values.
    ordered = np.sort(values)
    zero_count = sample_count - len(ordered)
    middle = (sample_count - 1) // 2, sample_count // 2
    median = (
        sum(
            0.0 if index < zero_count else float(ordered[index - zero_count])
            for index in middle
        )
        / 2
    )
    return {"total": total, "mean": total / sample_count, "median": median}


def _visualization_data(
    feature_load: biom.Table,
    gene_taxonomy: pd.DataFrame,
    metadata: Metadata | None,
) -> dict:
    sample_ids = list(feature_load.ids(axis="sample"))
    matrix = feature_load.matrix_data.tocsr()
    mapping = _align_gene_taxonomy(feature_load, gene_taxonomy)
    decoded_pairs = list(
        mapping[["Taxon ID", "Gene ID"]].itertuples(index=False, name=None)
    )
    labels = mapping.drop_duplicates("Taxon ID").set_index("Taxon ID")[["Taxon"]]
    full_taxa, short_taxa = _taxon_labels(decoded_pairs, labels)
    pairs = []
    loads = []

    for row, (taxon_id, function_id) in enumerate(decoded_pairs):
        start, end = matrix.indptr[row], matrix.indptr[row + 1]
        values = matrix.data[start:end]
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("TFA loads must be finite and nonnegative.")
        pairs.append(
            {
                "taxon": full_taxa[taxon_id],
                "taxon_short": short_taxa[taxon_id],
                "taxon_id": taxon_id,
                "lineage": [
                    part.strip() for part in labels.at[taxon_id, "Taxon"].split(";")
                ],
                "function": function_id,
                **_sample_summaries(values, len(sample_ids)),
            }
        )
        loads.extend(
            [row, int(sample), float(value)]
            for sample, value in zip(matrix.indices[start:end], values)
            if value != 0
        )

    groups = {}
    if metadata is not None:
        frame = metadata.to_dataframe()
        for column, properties in metadata.columns.items():
            if properties.type == "categorical":
                series = frame[column].reindex(sample_ids)
                groups[column] = [
                    "(missing)" if pd.isna(value) else str(value) for value in series
                ]

    return {"samples": sample_ids, "pairs": pairs, "loads": loads, "groups": groups}


def explore_tfa(
    output_dir: str,
    feature_load: biom.Table,
    gene_taxonomy: pd.DataFrame,
    metadata: Metadata | None = None,
) -> None:
    """Write a TFA explorer using mapping labels and optional sample groups."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    assets = resources.files("q2_annotate") / "assets" / "tfa_explore"
    payload = json.dumps(
        _visualization_data(feature_load, gene_taxonomy, metadata),
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    )
    # JSON inside a script element must not contain literal HTML delimiters.
    payload = (
        payload.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    )
    html = (assets / "index.html").read_text(encoding="utf-8")
    (output / "index.html").write_text(
        html.replace("__TFA_DATA__", payload), encoding="utf-8"
    )
    for filename in ("explore.js", "taxonomy.js", "style.css"):
        shutil.copyfile(assets / filename, output / filename)
