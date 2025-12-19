/**
 * BUSCO Detailed View - Unified JavaScript Module
 * Handles both sample_data and feature_data detailed views
 */

// Global state
window.currentSortCategory = 'default';
window.currentSortSource = 'default';
window.currentSortOrder = 'desc';
window.sampleViews = {};
window.selectedSamples = new Set();
window.sampleColors = {};

// Viridis color palette for samples
const viridisColors = [
  '#440154', '#482777', '#3f4a8a', '#31688e', '#26828e',
  '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725',
];

/**
 * Calculate sort order map for MAG IDs based on specified category and source
 */
function getSortOrderMap(sampleBuscoData, sampleAssemblyData, category, source, sortOrder) {
  const magIds = [...new Set(sampleBuscoData.map(d => d.mag_id))].sort();

  let sortedMagIds;
  if (category === 'default' || source === 'default') {
    sortedMagIds = magIds;
  } else if (source === 'busco') {
    const valueMap = {};
    sampleBuscoData.forEach(d => {
      if (d.category === category) {
        valueMap[d.mag_id] = d.BUSCO_percentage || 0;
      }
    });
    sortedMagIds = [...magIds].sort((a, b) => {
      const valA = valueMap[a] || 0;
      const valB = valueMap[b] || 0;
      return sortOrder === 'asc' ? (valA - valB) : (valB - valA);
    });
  } else if (source === 'assembly') {
    const valueMap = {};
    sampleAssemblyData.forEach(d => {
      valueMap[d.mag_id] = d[category] || 0;
    });
    sortedMagIds = [...magIds].sort((a, b) => {
      const valA = valueMap[a] || 0;
      const valB = valueMap[b] || 0;
      return sortOrder === 'asc' ? (valA - valB) : (valB - valA);
    });
  } else {
    sortedMagIds = magIds;
  }

  // Create a map of mag_id -> sort_order
  const orderMap = {};
  sortedMagIds.forEach((magId, index) => {
    orderMap[magId] = index;
  });
  return orderMap;
}

/**
 * Get sorted MAG IDs (for feature_data view)
 */
function getSortedMagIds(detailedData, assemblyData, magIdsDefault, category, source) {
  if (category === 'default' || source === 'default') {
    return [...magIdsDefault];
  }

  const valueMap = {};

  if (source === 'busco') {
    detailedData.forEach(d => {
      if (d.category === category) {
        valueMap[d.mag_id] = d.BUSCO_percentage || 0;
      }
    });
    return [...magIdsDefault].sort((a, b) => {
      const valA = valueMap[a] || 0;
      const valB = valueMap[b] || 0;
      return valB - valA;  // Descending order
    });
  } else if (source === 'assembly') {
    assemblyData.forEach(d => {
      valueMap[d.mag_id] = d[category] || 0;
    });
    const ascending = (category === 'percent_gaps' || category === 'scaffolds');
    return [...magIdsDefault].sort((a, b) => {
      const valA = valueMap[a] || 0;
      const valB = valueMap[b] || 0;
      return ascending ? (valA - valB) : (valB - valA);
    });
  }

  return [...magIdsDefault];
}

/**
 * Add sort_order field to data based on order map
 */
function addSortOrder(data, orderMap) {
  return data.map(d => ({
    ...d,
    sort_order: orderMap[d.mag_id] !== undefined ? orderMap[d.mag_id] : 9999
  }));
}

/**
 * Update sort button styling to reflect active state
 */
function updateSortButtons(activeCategory) {
  document.querySelectorAll('.sort-btn').forEach(btn => {
    if (btn.dataset.sort === activeCategory) {
      btn.classList.remove('btn-outline-secondary');
      btn.classList.add('btn-primary');
    } else {
      btn.classList.remove('btn-primary');
      btn.classList.add('btn-outline-secondary');
    }
  });
}

/**
 * Toggle sample selection (sample_data view)
 */
function toggleSample(sampleId) {
  const button = document.querySelector(`[data-sample-id="${sampleId}"]`);

  if (window.selectedSamples.has(sampleId)) {
    window.selectedSamples.delete(sampleId);
    button.classList.remove('selected');
  } else {
    window.selectedSamples.add(sampleId);
    button.classList.add('selected');
  }

  updateSelectionCount();
  updatePlotCards();
}

/**
 * Select all samples (sample_data view)
 */
function selectAllSamples() {
  const buttons = document.querySelectorAll('.sample-selection-item');
  buttons.forEach(button => {
    const sampleId = button.dataset.sampleId;
    window.selectedSamples.add(sampleId);
    button.classList.add('selected');
  });
  updateSelectionCount();
  updatePlotCards();
}

/**
 * Clear all sample selections (sample_data view)
 */
function clearAllSamples() {
  const buttons = document.querySelectorAll('.sample-selection-item');
  buttons.forEach(button => {
    const sampleId = button.dataset.sampleId;
    window.selectedSamples.delete(sampleId);
    button.classList.remove('selected');
  });
  window.selectedSamples.clear();
  updateSelectionCount();
  updatePlotCards();
}

/**
 * Update selection count display (sample_data view)
 */
function updateSelectionCount() {
  const count = window.selectedSamples.size;
  const countElement = document.getElementById('selectionCount');
  if (countElement) {
    countElement.textContent = `${count} sample${count !== 1 ? 's' : ''} selected`;
  }
}

/**
 * Update plot cards container (sample_data view)
 */
function updatePlotCards() {
  const container = document.getElementById('plotCardsContainer');
  if (!container) return;

  container.innerHTML = '';
  window.sampleViews = {};

  window.selectedSamples.forEach(sampleId => {
    const card = createSamplePlotCard(sampleId);
    container.appendChild(card);
  });
}

/**
 * Create a plot card for a single sample (sample_data view)
 */
function createSamplePlotCard(sampleId) {
  const color = window.sampleColors[sampleId];

  const metricSelector = document.getElementById('metricSelector');
  const selectedMetric = metricSelector ? metricSelector.value : 'contigs_n50';
  const metricTitle = window.assemblyMetrics[selectedMetric] || selectedMetric;

  const sampleBuscoData = window.detailedData.filter(d => d.sample_id === sampleId);
  const sampleAssemblyData = window.assemblyData.filter(d => d.sample_id === sampleId);

  const cardDiv = document.createElement('div');
  cardDiv.className = 'col-lg-6 mb-2';
  cardDiv.id = `card-${sampleId}`;
  cardDiv.innerHTML = `
    <div class="card h-100">
      <div class="card-header d-flex justify-content-between align-items-center">
        <h6 class="mb-0">
          <span class="sample-color-indicator me-2" style="background-color: ${color};"></span>
          ${sampleId}
        </h6>
        <button class="btn btn-sm btn-outline-danger" onclick="removeSample('${sampleId}')">
          <i class="fas fa-times"></i> Hide
        </button>
      </div>
      <div class="card-body">
        <div class="row">
          <div class="col-8 pe-0">
            <div id="busco-plot-${sampleId}" class="plot-container">Loading BUSCO plot...</div>
          </div>
          <div class="col-4 ps-0">
            <div id="assembly-plot-${sampleId}" class="plot-container">Loading assembly plot...</div>
          </div>
        </div>
      </div>
    </div>
  `;

  if (!window.sampleViews[sampleId]) {
    window.sampleViews[sampleId] = { busco: null, assembly: null };
  }

  if (sampleBuscoData.length > 0) {
    setTimeout(async () => {
      const buscoContainer = document.getElementById(`busco-plot-${sampleId}`);
      const assemblyContainer = document.getElementById(`assembly-plot-${sampleId}`);

      if (!buscoContainer || !assemblyContainer) return;

      const orderMap = getSortOrderMap(sampleBuscoData, sampleAssemblyData, window.currentSortCategory, window.currentSortSource, window.currentSortOrder);
      const buscoDataWithSort = addSortOrder(sampleBuscoData, orderMap);
      const assemblyDataWithSort = addSortOrder(sampleAssemblyData, orderMap);

      const buscoPlotSpec = JSON.parse(JSON.stringify(window.buscoSpec));
      buscoPlotSpec.params[0].value = sampleId;
      buscoPlotSpec.data = { name: 'source', values: buscoDataWithSort };

      const assemblyPlotSpec = JSON.parse(JSON.stringify(window.assemblySpec));
      assemblyPlotSpec.params[0].value = sampleId;
      assemblyPlotSpec.params[1].value = selectedMetric;
      assemblyPlotSpec.params[2].value = metricTitle;
      assemblyPlotSpec.data = { name: 'source', values: assemblyDataWithSort };
      if (assemblyPlotSpec.encoding && assemblyPlotSpec.encoding.x) {
        assemblyPlotSpec.encoding.x.title = metricTitle;
      }

      try {
        const buscoResult = await vegaEmbed(buscoContainer, buscoPlotSpec, {
          actions: false,
          renderer: 'svg'
        });
        buscoResult.view.logLevel(vega.Warn);
        window.sampleViews[sampleId].busco = buscoResult.view;
        window.sampleViews[sampleId].buscoData = sampleBuscoData;
      } catch (error) {
        buscoContainer.innerHTML = `<div class="alert alert-warning">Error: ${error.message}</div>`;
      }

      try {
        const assemblyResult = await vegaEmbed(assemblyContainer, assemblyPlotSpec, {
          actions: false,
          renderer: 'svg'
        });
        assemblyResult.view.logLevel(vega.Warn);
        window.sampleViews[sampleId].assembly = assemblyResult.view;
        window.sampleViews[sampleId].assemblyData = sampleAssemblyData;
      } catch (error) {
        assemblyContainer.innerHTML = `<div class="alert alert-warning">Error: ${error.message}</div>`;
      }
    }, 100);
  } else {
    setTimeout(() => {
      const buscoContainer = document.getElementById(`busco-plot-${sampleId}`);
      const assemblyContainer = document.getElementById(`assembly-plot-${sampleId}`);
      if (buscoContainer) {
        buscoContainer.innerHTML = '<div class="alert alert-info">No BUSCO data available.</div>';
      }
      if (assemblyContainer) {
        assemblyContainer.innerHTML = '<div class="alert alert-info">No assembly data available.</div>';
      }
    }, 100);
  }

  return cardDiv;
}

/**
 * Update assembly metric for all views (sample_data view)
 */
async function updateAssemblyMetric(newMetric) {
  const metricTitle = window.assemblyMetrics[newMetric] || newMetric;

  for (const sampleId of window.selectedSamples) {
    const views = window.sampleViews[sampleId];
    if (!views) continue;

    const sampleAssemblyData = views.assemblyData || [];
    if (sampleAssemblyData.length === 0) continue;

    const assemblyContainer = document.getElementById(`assembly-plot-${sampleId}`);
    if (!assemblyContainer) continue;

    // Calculate current sort order
    const sampleBuscoData = views.buscoData || [];
    const orderMap = getSortOrderMap(sampleBuscoData, sampleAssemblyData, window.currentSortCategory, window.currentSortSource, window.currentSortOrder);
    const assemblyDataWithSort = addSortOrder(sampleAssemblyData, orderMap);

    // Create new spec with updated metric
    const assemblyPlotSpec = JSON.parse(JSON.stringify(window.assemblySpec));
    assemblyPlotSpec.params[0].value = sampleId;
    assemblyPlotSpec.params[1].value = newMetric;
    assemblyPlotSpec.params[2].value = metricTitle;
    assemblyPlotSpec.data = { name: 'source', values: assemblyDataWithSort };
    if (assemblyPlotSpec.encoding && assemblyPlotSpec.encoding.x) {
      assemblyPlotSpec.encoding.x.title = metricTitle;
    }

    try {
      const assemblyResult = await vegaEmbed(assemblyContainer, assemblyPlotSpec, {
        actions: false,
        renderer: 'svg'
      });
      assemblyResult.view.logLevel(vega.Warn);
      views.assembly = assemblyResult.view;
    } catch (error) {
      console.error('Error updating assembly plot for', sampleId, ':', error);
    }
  }
}

/**
 * Update sort order for all views (sample_data view)
 */
function updateSortOrder() {
  window.selectedSamples.forEach(sampleId => {
    const views = window.sampleViews[sampleId];
    if (!views) return;

    const sampleBuscoData = views.buscoData || [];
    const sampleAssemblyData = views.assemblyData || [];

    if (sampleBuscoData.length === 0) return;

    const orderMap = getSortOrderMap(sampleBuscoData, sampleAssemblyData, window.currentSortCategory, window.currentSortSource, window.currentSortOrder);

    if (views.busco) {
      try {
        const buscoDataWithSort = addSortOrder(sampleBuscoData, orderMap);
        const cs = vega.changeset().remove(() => true).insert(buscoDataWithSort);
        views.busco.change('source', cs).run();
      } catch (e) {
        console.error('Error updating BUSCO view for', sampleId, ':', e);
      }
    }

    if (views.assembly) {
      try {
        const assemblyDataWithSort = addSortOrder(sampleAssemblyData, orderMap);
        const cs = vega.changeset().remove(() => true).insert(assemblyDataWithSort);
        views.assembly.change('source', cs).run();
      } catch (e) {
        console.error('Error updating Assembly view for', sampleId, ':', e);
      }
    }
  });
}

/**
 * Remove sample from view (sample_data view)
 */
function removeSample(sampleId) {
  const button = document.querySelector(`[data-sample-id="${sampleId}"]`);
  window.selectedSamples.delete(sampleId);
  button.classList.remove('btn-primary');
  button.classList.add('btn-outline-secondary');

  if (window.sampleViews[sampleId]) {
    delete window.sampleViews[sampleId];
  }

  const card = document.getElementById(`card-${sampleId}`);
  if (card) {
    card.remove();
  }

  updateSelectionCount();
}

/**
 * Update all plots (sample_data view)
 */
function updateAllPlots() {
  updatePlotCards();
}

/**
 * Update plots (feature_data view)
 */
function updatePlots() {
  const metricSelector = document.getElementById('metricSelector');
  const selectedMetric = metricSelector ? metricSelector.value : 'contigs_n50';
  const metricTitle = window.assemblyMetrics[selectedMetric] || selectedMetric;

  const buscoContainer = document.getElementById('buscoPlotContainer');
  const assemblyContainer = document.getElementById('assemblyPlotContainer');

  if (!buscoContainer || !assemblyContainer) return;

  let magIdsSorted = getSortedMagIds(window.detailedData, window.assemblyData, window.magIdsDefault, window.currentSortCategory, window.currentSortSource);

  const buscoPlotSpec = JSON.parse(JSON.stringify(window.buscoSpec));
  buscoPlotSpec.data = { values: window.detailedData };
  if (buscoPlotSpec.encoding && buscoPlotSpec.encoding.y) {
    buscoPlotSpec.encoding.y.sort = magIdsSorted;
  }

  const assemblyPlotSpec = JSON.parse(JSON.stringify(window.assemblySpec));
  assemblyPlotSpec.params[1].value = selectedMetric;
  assemblyPlotSpec.params[2].value = metricTitle;
  assemblyPlotSpec.data = { values: window.assemblyData };
  if (assemblyPlotSpec.encoding && assemblyPlotSpec.encoding.x) {
    assemblyPlotSpec.encoding.x.title = metricTitle;
  }
  if (assemblyPlotSpec.encoding && assemblyPlotSpec.encoding.y) {
    assemblyPlotSpec.encoding.y.sort = magIdsSorted;
  }

  vegaEmbed(buscoContainer, buscoPlotSpec, {
    actions: false,
    renderer: 'svg'
  }).then(result => {
    result.view.logLevel(vega.Warn);
  }).catch(error => {
    buscoContainer.innerHTML = `<div class="alert alert-warning">Error: ${error.message}</div>`;
  });

  vegaEmbed(assemblyContainer, assemblyPlotSpec, {
    actions: false,
    renderer: 'svg'
  }).then(result => {
    result.view.logLevel(vega.Warn);
  }).catch(error => {
    assemblyContainer.innerHTML = `<div class="alert alert-warning">Error: ${error.message}</div>`;
  });
}

/**
 * Initialize sample_data view
 */
function initSampleDataView() {
  const detailedData = JSON.parse(document.getElementById('detailed_data').textContent);
  const assemblyData = JSON.parse(document.getElementById('assembly_data').textContent);
  const sampleIds = JSON.parse(document.getElementById('sample_ids_json').textContent);
  const assemblyMetrics = JSON.parse(document.getElementById('assembly_metrics_json').textContent);

  const buscoSpec = JSON.parse(document.getElementById('vega_busco_detailed_spec').textContent);
  const assemblySpec = JSON.parse(document.getElementById('vega_assembly_detailed_spec').textContent);

  window.detailedData = detailedData;
  window.assemblyData = assemblyData;
  window.buscoSpec = buscoSpec;
  window.assemblySpec = assemblySpec;
  window.assemblyMetrics = assemblyMetrics;

  const sampleSelectionBox = document.getElementById('sampleSelectionBox');
  sampleIds.forEach((sampleId, index) => {
    const color = viridisColors[index % viridisColors.length];
    window.sampleColors[sampleId] = color;

    const sampleButton = document.createElement('div');
    const isFirstOrSecondSample = index === 0 || index === 1;
    sampleButton.className = `sample-selection-item ${isFirstOrSecondSample ? 'selected' : ''}`;
    sampleButton.innerHTML = `<span>${sampleId}</span>`;
    sampleButton.dataset.sampleId = sampleId;
    sampleButton.onclick = () => toggleSample(sampleId);

    if (isFirstOrSecondSample) {
      window.selectedSamples.add(sampleId);
    }

    sampleSelectionBox.appendChild(sampleButton);
  });

  document.getElementById('selectAllBtn').onclick = selectAllSamples;
  document.getElementById('clearAllBtn').onclick = clearAllSamples;

  const metricSelector = document.getElementById('metricSelector');
  if (metricSelector) {
    metricSelector.onchange = function() {
      updateAssemblyMetric(this.value);
    };
  }

  document.querySelectorAll('.sort-btn').forEach(btn => {
    btn.addEventListener('click', function() {
      window.currentSortCategory = this.dataset.sort;
      window.currentSortSource = this.dataset.source || 'default';
      updateSortButtons(window.currentSortCategory);

      // If sorting by assembly metric, update the metric selector and displayed metric
      if (window.currentSortSource === 'assembly' && metricSelector) {
        metricSelector.value = window.currentSortCategory;
        updateAssemblyMetric(window.currentSortCategory);
      }

      updateSortOrder();
    });
  });

  const sortOrderSelector = document.getElementById('sortOrderSelector');
  if (sortOrderSelector) {
    sortOrderSelector.addEventListener('change', function() {
      window.currentSortOrder = this.value;
      if (window.currentSortCategory !== 'default') {
        updateSortOrder();
      }
    });
  }

  updateSelectionCount();
  updatePlotCards();
}

/**
 * Initialize feature_data view
 */
function initFeatureDataView() {
  const buscoSpec = JSON.parse(document.getElementById('vega_busco_detailed_spec').textContent);
  const assemblySpec = JSON.parse(document.getElementById('vega_assembly_detailed_spec').textContent);
  const detailedData = JSON.parse(document.getElementById('detailed_data').textContent);
  const assemblyData = JSON.parse(document.getElementById('assembly_data').textContent);
  const magIdsDefault = JSON.parse(document.getElementById('mag_ids_sorted').textContent);
  const assemblyMetrics = JSON.parse(document.getElementById('assembly_metrics_json').textContent);

  window.buscoSpec = buscoSpec;
  window.assemblySpec = assemblySpec;
  window.detailedData = detailedData;
  window.assemblyData = assemblyData;
  window.magIdsDefault = magIdsDefault;
  window.assemblyMetrics = assemblyMetrics;

  const metricSelector = document.getElementById('metricSelector');
  if (metricSelector) {
    metricSelector.addEventListener('change', updatePlots);
  }

  document.querySelectorAll('.sort-btn').forEach(btn => {
    btn.addEventListener('click', function() {
      window.currentSortCategory = this.dataset.sort;
      window.currentSortSource = this.dataset.source || 'default';
      updateSortButtons(window.currentSortCategory);

      if (window.currentSortSource === 'assembly' && metricSelector) {
        metricSelector.value = window.currentSortCategory;
      }

      updatePlots();
    });
  });

  updatePlots();
}

/**
 * Initialize the appropriate view based on DOM elements present
 */
$(document).ready(function () {
  adjustTagsToBS3();

  // Detect which view to initialize based on DOM elements
  if (document.getElementById('sampleSelectionBox')) {
    initSampleDataView();
  } else if (document.getElementById('buscoPlotContainer')) {
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
