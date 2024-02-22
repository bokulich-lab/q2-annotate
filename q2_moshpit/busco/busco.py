# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import os
import tempfile
import q2_moshpit.busco.utils
from q2_moshpit.busco.utils import (
    _parse_busco_params,
    _render_html,
)
from q2_moshpit._utils import _process_common_input_params
from q2_types.per_sample_sequences._format import MultiMAGSequencesDirFmt


def evaluate_busco(
    output_dir: str,
    bins: MultiMAGSequencesDirFmt,
    mode: str = "genome",
    lineage_dataset: str = None,
    augustus: bool = False,
    augustus_parameters: str = None,
    augustus_species: str = None,
    auto_lineage: bool = False,
    auto_lineage_euk: bool = False,
    auto_lineage_prok: bool = False,
    cpu: int = 1,
    config: str = None,
    contig_break: int = 10,
    evalue: float = 1e-03,
    force: bool = False,
    limit: int = 3,
    long: bool = False,
    metaeuk_parameters: str = None,
    metaeuk_rerun_parameters: str = None,
    miniprot: bool = False,
    scaffold_composition: bool = False,
) -> None:
    """
    qiime2 visualization for the BUSCO assessment tool
    <https://busco.ezlab.org/>.

    Args:
        see all possible inputs by running `qiime moshpit plot_busco`

    Output:
        plots.zip: zip file containing all of the busco plots
        busco_output: all busco output files
        qiime_html: html for rendering the output plots
    """

    # Create dictionary with local variables
    # (kwargs passed to the function or their defaults) excluding
    # "output_dir" and "bins"
    kwargs = {
        k: v for k, v in locals().items() if k not in ["output_dir", "bins"]
    }

    # Filter out all kwargs that are None, False or 0.0
    common_args = _process_common_input_params(
        processing_func=_parse_busco_params, params=kwargs
    )

    # Creates output directory with path 'tmp'
    with tempfile.TemporaryDirectory() as tmp:
        # Run busco for every sample. Returns dictionary to report files.
        # Result NOT included in final output
        busco_results_dir = os.path.join(tmp, "busco_output")
        path_to_run_summaries = q2_moshpit.busco.utils._run_busco(
            output_dir=busco_results_dir,
            mags=bins,
            params=common_args,
        )

        # Collect result for each sample and save to file.
        # Result included in final output (file for download)
        all_summaries_path = os.path.join(
            output_dir, "all_batch_summaries.csv"
        )
        all_summaries_df = q2_moshpit.busco.utils._collect_summaries_and_save(
            all_summaries_path=all_summaries_path,
            path_to_run_summaries=path_to_run_summaries,
        )

        # Draw BUSCO plots for all samples
        # Result NOT included in final output
        plots_dir = os.path.join(tmp, "plots")
        paths_to_plots = q2_moshpit.busco.utils._draw_busco_plots(
            path_to_run_summaries=path_to_run_summaries,
            plots_dir=plots_dir
        )

        # Zip graphs for user download
        # Result included in final output (file for download)
        zip_name = os.path.join(output_dir, "busco_plots.zip")
        q2_moshpit.busco.utils._zip_busco_plots(
            paths_to_plots=paths_to_plots,
            zip_path=zip_name
        )

        # Render qiime html report
        # Result included in final output
        _render_html(output_dir, all_summaries_df)
