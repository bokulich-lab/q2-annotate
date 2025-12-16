/**
 * BUSCO Index View - Unified JavaScript Module
 * Handles both sample_data and feature_data index views
 */

// Global state
window.views = [];
window.currentViews = [];

// Configuration
const actionsCfg = { actions: false };

// Metric titles for display
const metricTitles = {
  'single': '% single',
  'duplicated': '% duplicated',
  'fragmented': '% fragmented',
  'missing': '% missing',
  'completeness': '% completeness',
  'contamination': '% contamination',
  'contigs_n50': 'contig N50 [bp]',
  'length': 'length [bp]'
};

// Y-axis formats by metric type
const yFormats = {
  'single': '.1f',
  'duplicated': '.1f',
  'fragmented': '.1f',
  'missing': '.1f',
  'completeness': '.1f',
  'contamination': '.1f',
  'contigs_n50': '.2s',
  'length': '.2s'
};

/**
 * Populate summary statistics table
 */
function populateSummaryTable(summarySpec, columns) {
  const rows = ['rowMinimum', 'rowMedian', 'rowMean', 'rowMaximum', 'rowTotal'];
  rows.forEach((rowId, index) => {
    let rowElement = document.getElementById(rowId);
    if (!rowElement) return;

    let colVals = summarySpec.data[index];
    columns.forEach((col) => {
      let td = document.createElement('td');
      td.textContent = colVals[col];
      rowElement.appendChild(td);
    });
  });
}

/**
 * Clear current plots
 */
function clearPlots(grid) {
  grid.innerHTML = '';
  window.currentViews.forEach(view => {
    const index = window.views.indexOf(view);
    if (index > -1) {
      window.views.splice(index, 1);
    }
  });
  window.currentViews = [];
}

/**
 * Render histograms
 */
function renderHistograms(grid, histogramSpec, histogramData, metrics) {
  clearPlots(grid);
  grid.className = 'square-chart-grid';

  metrics.forEach((metric) => {
    const div = document.createElement('div');
    grid.appendChild(div);

    const spec = JSON.parse(JSON.stringify(histogramSpec));
    spec.data = { values: histogramData };
    spec.params[1].value = metric;
    spec.encoding.x.title = metricTitles[metric] || metric;

    vegaEmbed(div, spec, actionsCfg).then(
      function (result) {
        result.view.logLevel(vega.Warn);
        result.view.resize();
        window.views.push(result.view);
        window.currentViews.push(result.view);
      }
    ).catch(
      function (error) {
        handleErrors([error], $(div));
      }
    );
  });
}

/**
 * Render box plots
 */
function renderBoxPlots(grid, boxPlotSpec, boxPlotData, metrics, hidePoints) {
  clearPlots(grid);
  grid.className = 'chart-grid';

  metrics.forEach((metric) => {
    if (!boxPlotData[metric]) return;

    const div = document.createElement('div');
    grid.appendChild(div);

    const spec = JSON.parse(JSON.stringify(boxPlotSpec));
    spec.data = { values: boxPlotData[metric] };
    spec.layer[0].encoding.y.axis.format = yFormats[metric] || '.2f';

    // Set hide_points param
    spec.params.forEach(param => {
      if (param.name === 'hide_points') {
        param.value = hidePoints;
      }
    });

    // Update y-axis title in layer encodings
    if (spec.layer) {
      spec.layer.forEach(layer => {
        if (layer.encoding && layer.encoding.y) {
          layer.encoding.y.title = metricTitles[metric] || metric;
        }
      });
    }

    vegaEmbed(div, spec, actionsCfg).then(
      function (result) {
        result.view.logLevel(vega.Warn);
        result.view.resize();
        window.views.push(result.view);
        window.currentViews.push(result.view);
      }
    ).catch(
      function (error) {
        handleErrors([error], $(div));
      }
    );
  });
}

/**
 * Update description text based on visualization type
 */
function updateDescription(descriptionEl, type, hasSampleFilter) {
  if (type === 'histograms') {
    if (hasSampleFilter) {
      descriptionEl.textContent = 'The histograms below show distributions of BUSCO marker count fractions and assembly metrics. Use the sample filter to show all samples or a single sample. Switch between histograms and box plots using the view toggle.';
    } else {
      descriptionEl.textContent = 'The histograms below show distributions of BUSCO marker count fractions and assembly metrics. Switch between histograms and box plots using the view toggle.';
    }
  } else {
    if (hasSampleFilter) {
      descriptionEl.textContent = 'The box plots below show distributions of BUSCO marker count fractions and assembly metrics. Each box shows the median, quartiles, and outliers for each category. Individual data points are overlaid with low transparency. Use the sample filter to highlight a specific sample across all box plots.';
    } else {
      descriptionEl.textContent = 'The box plots below show distributions of BUSCO marker count fractions and assembly metrics. Each box shows the median, quartiles, and outliers for each category. Individual data points are overlaid with low transparency.';
    }
  }
}

/**
 * Update all views with selected ID signal
 */
function updateViewsSelectedId(value) {
  window.views.forEach(v => {
    try { v.signal('selected_id', value).run(); } catch (e) {}
  });
}

/**
 * Embed scatter plot
 */
function embedScatterPlot(containerSelector, scatterSpec, scatterData, upperX, upperY, colorField, colorTitle, filterTransform) {
  scatterSpec.data = { values: scatterData };
  scatterSpec.params[1].value = upperX;
  scatterSpec.params[2].value = upperY;

  if (colorField) {
    scatterSpec.encoding.color.field = colorField;
    scatterSpec.encoding.color.title = colorTitle;
  }

  if (filterTransform) {
    scatterSpec.transform[0].filter = filterTransform;
  }

  vegaEmbed(containerSelector, scatterSpec, actionsCfg).then(
    function (result) {
      result.view.logLevel(vega.Warn);
      result.view.resize();
      window.views.push(result.view);
    }
  ).catch(
    function (error) {
      handleErrors([error], $(containerSelector));
    }
  );
}

/**
 * Embed unbinned plot (sample_data only)
 */
function embedUnbinnedPlot(containerSelector, unbinnedSpec, unbinnedData) {
  unbinnedSpec = JSON.parse(JSON.stringify(unbinnedSpec));
  unbinnedSpec.data = {
    values: unbinnedData.map(d => ({
      ...d,
      metric: d.unbinned_contigs_count,
      category: 'unbinned_contigs_count'
    }))
  };
  unbinnedSpec.params[1].value = 'unbinned_contigs_count';
  unbinnedSpec.encoding.x.title = 'Unbinned contig count';
  unbinnedSpec.encoding.y.title = 'Sample count';

  vegaEmbed(containerSelector, unbinnedSpec, actionsCfg).then(
    function (result) {
      result.view.logLevel(vega.Warn);
      result.view.resize();
    }
  ).catch(
    function (error) {
      handleErrors([error], $(containerSelector));
    }
  );
}

/**
 * Initialize sample_data view
 */
function initSampleDataView() {
  const summarySpec = JSON.parse(document.getElementById('summary_stats_json').textContent);
  const columns = ['single', 'duplicated', 'fragmented', 'missing', 'complete', 'completeness', 'contamination', 'unbinned_contigs'];
  populateSummaryTable(summarySpec, columns);

  // Populate sample dropdown
  const selectEl = document.getElementById('globalSampleSelect');
  try {
    const ids = JSON.parse(document.getElementById('sample_ids_json').textContent || '[]');
    ids.forEach(id => {
      const opt = document.createElement('option');
      opt.value = id;
      opt.textContent = id;
      selectEl.appendChild(opt);
    });
  } catch (e) {}

  // Load specs and data
  const histogramSpec = JSON.parse(document.getElementById('vega_histogram_spec').textContent);
  const boxPlotSpec = JSON.parse(document.getElementById('vega_box_plot_spec').textContent);
  const histogramData = JSON.parse(document.getElementById('histogram_data').textContent);
  const boxPlotData = JSON.parse(document.getElementById('box_plot_data').textContent);
  const metrics = JSON.parse(document.getElementById('metrics_json').textContent);

  const grid = document.getElementById('plotMarkers');
  const visualizationToggle = document.getElementById('visualizationToggle');
  const descriptionEl = document.getElementById('visualizationDescription');
  const hidePointsContainer = document.getElementById('hidePointsContainer');
  const hidePointsCheckbox = document.getElementById('hidePointsCheckbox');

  // Toggle event listener
  visualizationToggle.addEventListener('change', function() {
    const selectedType = this.value;
    updateDescription(descriptionEl, selectedType, true);
    hidePointsContainer.style.display = (selectedType === 'boxplots') ? 'flex' : 'none';

    if (selectedType === 'histograms') {
      renderHistograms(grid, histogramSpec, histogramData, metrics);
    } else {
      renderBoxPlots(grid, boxPlotSpec, boxPlotData, metrics, hidePointsCheckbox.checked);
    }
  });

  // Hide points checkbox listener
  hidePointsCheckbox.addEventListener('change', function() {
    if (visualizationToggle.value === 'boxplots') {
      renderBoxPlots(grid, boxPlotSpec, boxPlotData, metrics, hidePointsCheckbox.checked);
    }
  });

  // Initial render (box plots by default)
  renderBoxPlots(grid, boxPlotSpec, boxPlotData, metrics, hidePointsCheckbox.checked);

  // Embed scatter plot if comp_cont is enabled
  const scatterSpecEl = document.getElementById('vega_scatter_spec');
  const scatterDataEl = document.getElementById('scatter_data');
  if (scatterSpecEl && scatterDataEl) {
    const scatterSpec = JSON.parse(scatterSpecEl.textContent);
    const scatterData = JSON.parse(scatterDataEl.textContent);
    const upperX = window.buscoUpperX || 100;
    const upperY = window.buscoUpperY || 100;
    embedScatterPlot("#scatterPlotDiv", scatterSpec, scatterData, upperX, upperY);
  }

  // Embed unbinned plot if unbinned data is available
  const unbinnedSpecEl = document.getElementById('vega_unbinned_spec');
  const unbinnedDataEl = document.getElementById('unbinned_data');
  if (unbinnedSpecEl && unbinnedDataEl) {
    const unbinnedSpec = JSON.parse(unbinnedSpecEl.textContent);
    const unbinnedData = JSON.parse(unbinnedDataEl.textContent);
    embedUnbinnedPlot("#plotUnbinned", unbinnedSpec, unbinnedData);
  }

  // Sample selector listener
  selectEl.addEventListener('change', (e) => {
    const value = e.target.value || 'All';
    updateViewsSelectedId(value);
  });

  // Initialize with 'All'
  updateViewsSelectedId('All');

  // Handle collapse chevron rotation
  const collapseElement = document.getElementById('histogramsCollapse');
  const chevronIcon = document.querySelector('[data-bs-target="#histogramsCollapse"] i');

  if (collapseElement && chevronIcon) {
    collapseElement.addEventListener('show.bs.collapse', function () {
      chevronIcon.classList.remove('fa-chevron-right');
      chevronIcon.classList.add('fa-chevron-down');
    });

    collapseElement.addEventListener('hide.bs.collapse', function () {
      chevronIcon.classList.remove('fa-chevron-down');
      chevronIcon.classList.add('fa-chevron-right');
    });
  }
}

/**
 * Initialize feature_data view
 */
function initFeatureDataView() {
  const summarySpec = JSON.parse(document.getElementById('summary_stats_json').textContent);
  const columns = ['single', 'duplicated', 'fragmented', 'missing', 'complete', 'completeness', 'contamination'];
  populateSummaryTable(summarySpec, columns);

  // Load specs and data
  const histogramSpec = JSON.parse(document.getElementById('vega_histogram_spec').textContent);
  const boxPlotSpec = JSON.parse(document.getElementById('vega_box_plot_spec').textContent);
  const histogramData = JSON.parse(document.getElementById('histogram_data').textContent);
  const boxPlotData = JSON.parse(document.getElementById('box_plot_data').textContent);
  const metrics = JSON.parse(document.getElementById('metrics_json').textContent);

  const grid = document.getElementById('plotMarkers');
  const visualizationToggle = document.getElementById('visualizationToggle');
  const descriptionEl = document.getElementById('visualizationDescription');
  const hidePointsContainer = document.getElementById('hidePointsContainer');
  const hidePointsCheckbox = document.getElementById('hidePointsCheckbox');

  // Toggle event listener
  visualizationToggle.addEventListener('change', function() {
    const selectedType = this.value;
    updateDescription(descriptionEl, selectedType, false);
    hidePointsContainer.style.display = (selectedType === 'boxplots') ? 'flex' : 'none';

    if (selectedType === 'histograms') {
      renderHistograms(grid, histogramSpec, histogramData, metrics);
    } else {
      renderBoxPlots(grid, boxPlotSpec, boxPlotData, metrics, hidePointsCheckbox.checked);
    }
  });

  // Hide points checkbox listener
  hidePointsCheckbox.addEventListener('change', function() {
    if (visualizationToggle.value === 'boxplots') {
      renderBoxPlots(grid, boxPlotSpec, boxPlotData, metrics, hidePointsCheckbox.checked);
    }
  });

  // Initial render (box plots by default)
  renderBoxPlots(grid, boxPlotSpec, boxPlotData, metrics, hidePointsCheckbox.checked);

  // Embed scatter plot if comp_cont is enabled
  const scatterSpecEl = document.getElementById('vega_scatter_spec');
  const scatterDataEl = document.getElementById('scatter_data');
  if (scatterSpecEl && scatterDataEl) {
    const scatterSpec = JSON.parse(scatterSpecEl.textContent);
    const scatterData = JSON.parse(scatterDataEl.textContent);
    const upperX = window.buscoUpperX || 100;
    const upperY = window.buscoUpperY || 100;
    embedScatterPlot(
      "#scatterPlotDiv",
      scatterSpec,
      scatterData,
      upperX,
      upperY,
      'mag_id',
      'MAG ID',
      "(selected_id == 'All') || (datum.mag_id == selected_id)"
    );
  }
}

/**
 * Initialize the appropriate view based on DOM elements present
 */
$(document).ready(function () {
  adjustTagsToBS3();

  // Detect which view to initialize based on presence of sample selector
  if (document.getElementById('globalSampleSelect')) {
    initSampleDataView();
  } else {
    initFeatureDataView();
  }
});

/**
 * Error handler for Vega/Vega-Lite
 */
function handleErrors(errors, element) {
  console.error('Vega/Vega-Lite error:', errors);
  if (element && element.length) {
    element.html('<div class="alert alert-danger"><strong>Error:</strong> Failed to render plot. Check console for details.</div>');
  }
}
