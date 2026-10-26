"use strict";

// Browser controls and Vega plots; taxonomy aggregation lives in taxonomy.js.
// data.json stores nonzero loads as [pair index, sample index, load] triplets.
let tfa;
const showUnclassified = document.getElementById("show-unclassified");
const levelSelect = document.getElementById("level");
const taxonSearch = document.getElementById("taxon-search");
const selectedTaxonomy = document.getElementById("selected-taxonomy");
const taxonSelect = document.getElementById("taxon");
const functionSelect = document.getElementById("function");
const groupSelect = document.getElementById("group");
const metricSelect = document.getElementById("metric");
const limitSelect = document.getElementById("limit");
const clearSelections = document.getElementById("clear-selections");
const status = document.getElementById("status");
const schema = "https://vega.github.io/schema/vega-lite/v5.json";
let activeViews = [];
let renderQueue = Promise.resolve();

function addOptions(select, labels, firstLabel) {
  const first = document.createElement("option");
  first.value = "";
  first.textContent = firstLabel;
  select.appendChild(first);
  for (const label of labels) {
    const option = document.createElement("option");
    option.value = label;
    option.textContent = label;
    select.appendChild(option);
  }
}

// Keep the original data in tfa and the selected taxonomy level in current.
let current;
let loadsByPair = [];
let deepestLevel = 1;
const controls = [showUnclassified, levelSelect, taxonSearch, taxonSelect,
                  functionSelect, groupSelect, metricSelect, limitSelect, clearSelections];
for (const control of controls) control.disabled = true;

/** Load data and populate controls before drawing the first plots. */
async function initializeExplorer() {
  const response = await fetch("data.json");
  if (!response.ok) throw new Error(`HTTP ${response.status}`);
  tfa = await response.json();
  const functions = [...new Set(tfa.pairs.map(pair => pair.function))].sort();
  deepestLevel = tfa.pairs.reduce((depth, pair) => Math.max(depth, pair.lineage.length), 1);
  for (let level = 1; level <= deepestLevel; level++) {
    levelSelect.add(new Option(`Level ${level}`, String(level)));
  }
  addOptions(functionSelect, functions, "All functions");
  addOptions(groupSelect, Object.keys(tfa.groups), "All samples");
  for (const control of controls) control.disabled = false;
  resetSelections();
  await renderQueue;
}

/** Restore the same control values and plots used on first load. */
function resetSelections() {
  levelSelect.value = String(deepestLevel);
  taxonSearch.value = "";
  functionSelect.value = "";
  groupSelect.value = "";
  metricSelect.value = "sum";
  limitSelect.value = "100";
  showUnclassified.checked = false;
  changeLevel();
  update();
}

function pairVisible(pair) {
  return showUnclassified.checked || !isUnclassifiedTaxon(pair, Number(levelSelect.value));
}

function pairLabel(pair) {
  return pair.taxon_label || latestTaxonLabel(pair.lineage || [], pair.taxon_id);
}

/** Rebuild the taxon list from search and function filters, keeping valid selections. */
function refreshTaxa() {
  const previous = taxonSelect.value;
  const query = taxonSearch.value.trim().toLocaleLowerCase();
  const taxa = [...new Set(current.pairs
    .filter(pairVisible)
    .filter(pair => !functionSelect.value || pair.function === functionSelect.value)
    .map(pair => pair.taxon))]
    .filter(taxon => taxon.toLocaleLowerCase().includes(query)).sort();
  taxonSelect.replaceChildren();
  addOptions(taxonSelect, taxa, query ? "Select a matching taxon" : "All taxa");
  const labelByTaxon = new Map(current.pairs.map(pair => [pair.taxon, pairLabel(pair)]));
  const terminal = taxon => labelByTaxon.get(taxon);
  const counts = new Map();
  for (const taxon of taxa) counts.set(terminal(taxon), (counts.get(terminal(taxon)) || 0) + 1);
  for (const option of taxonSelect.options) {
    if (!option.value) continue;
    const label = terminal(option.value);
    const parent = option.value.split(";").slice(0, -1).join("; ").trim();
    option.textContent = counts.get(label) > 1 ? `${label} [${parent}]` : label;
    option.title = option.value;
  }
  taxonSelect.value = taxa.includes(previous) ? previous : "";
  showSelectedTaxonomy();
}

function showSelectedTaxonomy() {
  taxonSelect.title = taxonSelect.value || "All taxa";
  selectedTaxonomy.textContent = taxonSelect.value
    ? `Selected lineage: ${taxonSelect.value}`
    : taxonSearch.value.trim() ? "Search filters the taxon list. Select a taxon to filter the plots." : "";
}

/** Collapse to the chosen level and index sparse loads for repeated filtering. */
function changeLevel() {
  current = collapseTaxonomy(tfa, Number(levelSelect.value));
  loadsByPair = current.pairs.map(() => []);
  for (const load of current.loads) loadsByPair[load[0]].push(load);
  taxonSelect.value = "";
  refreshTaxa();
}

/** Find pair indices matching the current taxon, function, and visibility choices. */
function matchedPairs() {
  const indices = [];
  current.pairs.forEach((pair, index) => {
    if (pairVisible(pair) && (!taxonSelect.value || pair.taxon === taxonSelect.value) &&
        (!functionSelect.value || pair.function === functionSelect.value)) {
      indices.push(index);
    }
  });
  return indices;
}

/** Sum selected pairs per sample, retaining zero loads and metadata groups. */
function sampleValues(pairIndices) {
  const values = new Float64Array(tfa.samples.length);
  for (const pairIndex of pairIndices) {
    for (const [, sampleIndex, load] of loadsByPair[pairIndex]) {
      values[sampleIndex] += load;
    }
  }
  const grouping = tfa.groups[groupSelect.value];
  return tfa.samples.map((sample, index) => ({
    sample,
    group: grouping ? grouping[index] : "All samples",
    load: values[index]
  }));
}

function plotWidth(id) {
  // Vega's fit-x sizing includes axes, legends, and padding in this width.
  return Math.max(1, Math.floor(document.getElementById(id).clientWidth));
}

/** Rank links by the chosen metric and build a Vega-Lite heatmap specification. */
function heatmapSpec(indices) {
  const field = {sum: "total", mean: "mean", median: "median"}[metricSelect.value];
  const label = {sum: "Sum", mean: "Mean", median: "Median"}[metricSelect.value];
  const ranked = [...indices].sort((a, b) => current.pairs[b][field] - current.pairs[a][field]);
  const topIndices = new Set(ranked.slice(0, Number(limitSelect.value)));
  const top = [...topIndices].map(index => current.pairs[index]);
  const taxaOrder = [...new Set(top.map(pair => pair.taxon))];
  const functionOrder = [...new Set(top.map(pair => pair.function))];
  // Full lineages identify rows; short labels must not merge taxa with identical names.
  const axisLabels = Object.fromEntries(top.map(pair => [pair.taxon, pairLabel(pair)]));
  const taxa = new Set(taxaOrder);
  const functions = new Set(functionOrder);
  const seen = new Map(taxaOrder.map(taxon => [taxon, new Set()]));
  const cells = [];
  for (const index of ranked) {
    const pair = current.pairs[index];
    if (!taxa.has(pair.taxon) || !functions.has(pair.function)) continue;
    seen.get(pair.taxon).add(pair.function);
    cells.push({
      ...pair,
      present: true,
      top_link: topIndices.has(index),
      state: pair[field] === 0 ? "Zero summary" : "Positive summary",
      nonzero_samples: (loadsByPair[index] || []).filter(load => load[2] > 0).length
    });
  }
  for (const taxon of taxaOrder) {
    for (const func of functionOrder) {
      if (!seen.get(taxon).has(func)) {
        cells.push({taxon, function: func, present: false, top_link: false,
                    state: "No estimated load"});
      }
    }
  }
  return {
    $schema: schema,
    data: {values: cells},
    width: plotWidth("heatmap"),
    autosize: {type: "fit-x", contains: "padding", resize: true},
    height: {step: 29},
    // Keep absent pairs visible even though their numeric summaries are missing.
    mark: {
      type: "rect", cursor: "pointer", invalid: null,
      tooltip: {
        expr: `{
          'Taxon': datum.taxon,
          'Function': datum.function,
          ${JSON.stringify(`${label} load`)}:
            datum.present ? datum.${field} : 'No estimated load'
        }`
      }
    },
    encoding: {
      x: {
        field: "function", type: "nominal", sort: functionOrder, title: "Function",
        axis: {labels: functionOrder.length <= 30}
      },
      y: {
        field: "taxon", type: "nominal", sort: taxaOrder, title: "Taxon",
        axis: {labelLimit: Math.floor(plotWidth("heatmap") * 0.3), labelExpr: `${JSON.stringify(axisLabels)}[datum.label]`}
      },
      color: {
        condition: [
          {test: "!datum.present", value: "#f3f3f3"},
          {test: `datum.${field} === 0`, value: "#ffffff"}
        ],
        field, type: "quantitative", title: `${label} load`,
        scale: {scheme: "blues", zero: true}
      },
      stroke: {condition: {test: "datum.top_link", value: "#800000"}, value: "#bdbdbd"},
      strokeWidth: {condition: {test: "datum.top_link", value: 1}, value: 0.5},
    },
    config: {view: {stroke: null}, legend: {orient: "bottom"}}
  };
}

function sampleGroupColor(samples) {
  // Use every sample, including categories absent from the top 100 bars.
  return {
    field: "group", type: "nominal", legend: null,
    scale: {domain: [...new Set(samples.map(sample => sample.group))].sort(), scheme: "tableau10"}
  };
}

function boxplotSpec(samples) {
  return {
    $schema: schema,
    data: {values: samples},
    width: plotWidth("boxplot"),
    autosize: {type: "fit-x", contains: "padding", resize: true},
    height: 300,
    mark: {type: "boxplot", extent: 1.5},
    encoding: {
      x: {field: "group", type: "nominal", sort: sampleGroupColor(samples).scale.domain, title: groupSelect.value || "Samples"},
      y: {field: "load", type: "quantitative", title: "TFA load", scale: {zero: true}},
      color: sampleGroupColor(samples)
    },
    config: {view: {stroke: null}}
  };
}

function sampleSpec(samples) {
  const top = [...samples].sort((a, b) => b.load - a.load).slice(0, 100);
  return {
    $schema: schema,
    data: {values: top},
    width: plotWidth("samples"),
    autosize: {type: "fit-x", contains: "padding", resize: true},
    height: 300,
    mark: {type: "bar"},
    encoding: {
      x: {field: "sample", type: "nominal", sort: "-y", title: "Sample", axis: {labels: top.length <= 30}},
      y: {field: "load", type: "quantitative", title: "TFA load", scale: {zero: true}},
      color: sampleGroupColor(samples),
      tooltip: [
        {field: "sample", type: "nominal", title: "Sample"},
        {field: "group", type: "nominal", title: groupSelect.value || "Group"},
        {field: "load", type: "quantitative", title: "TFA load", format: ".4g"}
      ]
    },
    config: {view: {stroke: null}}
  };
}

/** Dispose of old Vega views and redraw all three plots for the current selection. */
async function draw() {
  for (const view of activeViews) view.finalize();
  activeViews = [];
  const indices = matchedPairs();
  const samples = sampleValues(indices);
  const method = metricSelect.value;
  status.textContent = `${indices.length} matching taxon–function pairs · ${samples.length} samples`;
  document.getElementById("heatmap-caption").textContent =
    `Each square is one pair at the selected taxonomy level; color shows its ${{sum: "sum", mean: "mean", median: "median"}[method]} load across all samples, including zeros. The top links by this metric select the axes; all intersections are shown. Dark borders mark the top links, white cells have a zero summary, and gray cells have no estimated load in the TFA table. Taxa are summed within each sample before this summary. Click a square to inspect that pair.`;
  document.getElementById("group-caption").textContent = groupSelect.value
    ? `Each box compares sample loads for one ${groupSelect.value} category. Zeros are included.`
    : "One box summarizes all samples. Add categorical metadata to compare groups.";

  if (!indices.length || !samples.length) {
    for (const id of ["heatmap-plot", "boxplot", "samples"]) {
      document.getElementById(id).textContent = "No values to plot for this selection.";
    }
    return;
  }
  if (typeof vegaEmbed !== "function") {
    status.textContent += " · Vega could not load; an internet connection is needed to display the plots.";
    return;
  }
  try {
    const heatmapContainer = document.getElementById("heatmap");
    heatmapContainer.scrollTop = 0;
    heatmapContainer.scrollLeft = 0;
    // Use fresh targets so a previous Vega view cannot hide the next render.
    for (const id of ["heatmap-plot", "boxplot", "samples"]) {
      document.getElementById(id).replaceChildren();
    }
    const options = {actions: false, renderer: "svg"};
    const [heatmap, boxplot, samplePlot] = await Promise.all([
      vegaEmbed("#heatmap-plot", heatmapSpec(indices), options),
      vegaEmbed("#boxplot", boxplotSpec(samples), options),
      vegaEmbed("#samples", sampleSpec(samples), options)
    ]);
    activeViews = [heatmap.view, boxplot.view, samplePlot.view];
    heatmap.view.addEventListener("click", (_, item) => {
      const pair = item && item.datum;
      if (pair && pair.present && pair.taxon && pair.function) {
        taxonSearch.value = "";
        functionSelect.value = pair.function;
        refreshTaxa();
        taxonSelect.value = pair.taxon;
        update();
      }
    });
  } catch (error) {
    status.textContent = `Could not draw the plots: ${error.message}`;
  }
}

function update() {
  if (!tfa) return;
  showSelectedTaxonomy();
  // Serialize asynchronous redraws so an older render cannot overwrite a newer one.
  renderQueue = renderQueue.then(draw).catch(error => {
    status.textContent = `Could not draw the plots: ${error.message}`;
  });
}

clearSelections.addEventListener("click", resetSelections);
// Rebuild the taxon choices whenever a control changes which taxa are available.
showUnclassified.addEventListener("change", () => {
  refreshTaxa();
  update();
});
levelSelect.addEventListener("change", () => {
  changeLevel();
  update();
});
taxonSearch.addEventListener("input", () => {
  refreshTaxa();
  update();
});
functionSelect.addEventListener("change", () => {
  refreshTaxa();
  update();
});
for (const select of [taxonSelect, groupSelect, metricSelect, limitSelect]) {
  select.addEventListener("change", update);
}
let resizeTimer;
window.addEventListener("resize", () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(update, 150);
});
const initialization = initializeExplorer().catch(error => {
  status.textContent = `Could not load TFA data: ${error.message}`;
});
