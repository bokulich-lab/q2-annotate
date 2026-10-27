"use strict";

// Browser controls and rendering; Vega-Lite templates are in plots.json.
// Taxonomy aggregation lives in taxonomy.js.
// data.json stores nonzero loads as [pair index, sample index, load] triplets.
let tfa;
let plotTemplates;
const showUnclassified = document.getElementById("show-unclassified");
const levelSelect = document.getElementById("level");
const taxonSearch = document.getElementById("taxon-search");
const selectedTaxonomy = document.getElementById("selected-taxonomy");
const taxonSelect = document.getElementById("taxon");
const functionSearch = document.getElementById("function-search");
const selectedFunction = document.getElementById("selected-function");
const functionSelect = document.getElementById("function");
const groupSelect = document.getElementById("group");
const metricSelect = document.getElementById("metric");
const limitSelect = document.getElementById("limit");
const clearSelections = document.getElementById("clear-selections");
const status = document.getElementById("status");
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
                  functionSearch, functionSelect, groupSelect, metricSelect, limitSelect, clearSelections];
for (const control of controls) control.disabled = true;

/** Load data and populate controls before drawing the first plots. */
async function initializeExplorer() {
  async function readJson(filename) {
    const response = await fetch(filename);
    if (!response.ok) throw new Error(`${filename}: HTTP ${response.status}`);
    return response.json();
  }
  [tfa, plotTemplates] = await Promise.all([
    readJson("data.json"), readJson("plots.json")
  ]);
  deepestLevel = tfa.pairs.reduce((depth, pair) => Math.max(depth, pair.lineage.length), 1);
  for (let level = 1; level <= deepestLevel; level++) {
    levelSelect.add(new Option(`Level ${level}`, String(level)));
  }
  addOptions(groupSelect, Object.keys(tfa.groups), "All samples");
  for (const control of controls) control.disabled = false;
  resetSelections();
  await renderQueue;
}

/** Restore the same control values and plots used on first load. */
function resetSelections() {
  levelSelect.value = String(deepestLevel);
  taxonSearch.value = "";
  functionSearch.value = "";
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

function taxonMatchesSearch(taxon) {
  return taxon.toLocaleLowerCase().includes(taxonSearch.value.trim().toLocaleLowerCase());
}

function functionMatchesSearch(func) {
  return func.toLocaleLowerCase().includes(functionSearch.value.trim().toLocaleLowerCase());
}

/** Rebuild the function list from visible heatmap links. */
function refreshFunctions(pairs) {
  const previous = functionSelect.value;
  const functions = [...new Set(pairs.map(pair => pair.function))].sort();
  functionSelect.replaceChildren();
  addOptions(functionSelect, functions, functionSearch.value.trim() ? "All matching functions" : "All functions");
  functionSelect.value = functions.includes(previous) ? previous : "";
}

/** Rebuild the taxon list from visible heatmap links. */
function refreshTaxa(pairs) {
  const previous = taxonSelect.value;
  const query = taxonSearch.value.trim().toLocaleLowerCase();
  const taxa = [...new Set(pairs.map(pair => pair.taxon))].sort();
  taxonSelect.replaceChildren();
  addOptions(taxonSelect, taxa, query ? "All matching taxa" : "All taxa");
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

/** Use the same top links for dropdown domains and heatmap axes. */
function refreshSelectors() {
  // Clearing an invalid selection can expose other matching links. Recompute
  // until both selections are stable; selections are only retained or cleared.
  while (true) {
    const previousTaxon = taxonSelect.value;
    const previousFunction = functionSelect.value;
    refreshLimitOptions();
    const {ranked} = heatmapLinks(matchedPairs());
    const pairs = ranked.map(index => current.pairs[index]);
    refreshFunctions(pairs);
    refreshTaxa(pairs);
    if (taxonSelect.value === previousTaxon && functionSelect.value === previousFunction) break;
  }
}

function showSelectedTaxonomy() {
  taxonSelect.title = taxonSelect.value || (taxonSearch.value.trim() ? "All matching taxa" : "All taxa");
  selectedTaxonomy.textContent = taxonSelect.value
    ? `Selected lineage: ${taxonSelect.value}`
    : taxonSearch.value.trim() ? "Search filters the taxon list and plots. Select a taxon to narrow the results further." : "";
}

function showSelectedFunction() {
  functionSelect.title = functionSelect.value || (functionSearch.value.trim() ? "All matching functions" : "All functions");
  selectedFunction.textContent = functionSelect.value
    ? `Selected function: ${functionSelect.value}`
    : functionSearch.value.trim() ? "Search filters the function list and plots. Select a function to narrow the results further." : "";
}

/** Collapse to the chosen level and index sparse loads for repeated filtering. */
function changeLevel() {
  current = collapseTaxonomy(tfa, Number(levelSelect.value));
  loadsByPair = current.pairs.map(() => []);
  for (const load of current.loads) loadsByPair[load[0]].push(load);
  taxonSelect.value = "";
}

/** Find pair indices matching the search, taxon, function, and visibility choices. */
function matchedPairs(ignoreSelections = false) {
  const indices = [];
  current.pairs.forEach((pair, index) => {
    if (pairVisible(pair) && taxonMatchesSearch(pair.taxon) && functionMatchesSearch(pair.function) &&
        (ignoreSelections || !taxonSelect.value || pair.taxon === taxonSelect.value) &&
        (ignoreSelections || !functionSelect.value || pair.function === functionSelect.value)) {
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

/** Rank matching links once using the metric and limit selected in the UI. */
function rankedLinks(indices) {
  const field = {sum: "total", mean: "mean", median: "median"}[metricSelect.value];
  const ranked = [...indices].sort((a, b) => current.pairs[b][field] - current.pairs[a][field]);
  const topIndices = new Set(ranked.slice(0, Number(limitSelect.value)));
  return {ranked, topIndices};
}

/** Disable cutoffs that would remove a selected taxon or function from the axes. */
function refreshLimitOptions() {
  const {ranked} = rankedLinks(matchedPairs(true));
  const firstTaxonRank = new Map(), firstFunctionRank = new Map();
  ranked.forEach((index, rank) => {
    const pair = current.pairs[index];
    if (!firstTaxonRank.has(pair.taxon)) firstTaxonRank.set(pair.taxon, rank + 1);
    if (!firstFunctionRank.has(pair.function)) firstFunctionRank.set(pair.function, rank + 1);
  });
  // Searches or taxonomy changes can remove a selection from the context entirely.
  if (taxonSelect.value && !firstTaxonRank.has(taxonSelect.value)) taxonSelect.value = "";
  if (functionSelect.value && !firstFunctionRank.has(functionSelect.value)) functionSelect.value = "";
  const requiredRank = Math.max(
    firstTaxonRank.get(taxonSelect.value) || 0,
    firstFunctionRank.get(functionSelect.value) || 0
  );
  for (const option of limitSelect.options) {
    option.disabled = Number(option.value) < requiredRank;
    option.title = option.disabled ? `The current selection requires at least the top ${requiredRank} links.` : "";
  }
  const selected = [...limitSelect.options].find(option => option.value === limitSelect.value);
  if (selected && selected.disabled) {
    const next = [...limitSelect.options].find(option => !option.disabled);
    if (next) limitSelect.value = next.value;
  }
}

/** Keep the reference ranking fixed while dropdown selections narrow its axes. */
function heatmapLinks(indices) {
  const {ranked, topIndices} = rankedLinks(matchedPairs(true));
  const top = [...topIndices].map(index => current.pairs[index]);
  const taxaOrder = [...new Set(top.map(pair => pair.taxon))];
  const functionOrder = [...new Set(top.map(pair => pair.function))];
  const taxa = new Set(taxaOrder), functions = new Set(functionOrder);
  const matched = new Set(indices);
  const visible = ranked.filter(index => {
    const pair = current.pairs[index];
    return matched.has(index) && taxa.has(pair.taxon) && functions.has(pair.function);
  });
  const visibleTaxa = new Set(visible.map(index => current.pairs[index].taxon));
  const visibleFunctions = new Set(visible.map(index => current.pairs[index].function));
  return {
    ranked: visible, topIndices,
    taxaOrder: taxaOrder.filter(taxon => visibleTaxa.has(taxon)),
    functionOrder: functionOrder.filter(func => visibleFunctions.has(func))
  };
}

/** Rank links by the chosen metric and build a Vega-Lite heatmap specification. */
function heatmapSpec(indices) {
  const field = {sum: "total", mean: "mean", median: "median"}[metricSelect.value];
  const label = {sum: "Sum", mean: "Mean", median: "Median"}[metricSelect.value];
  const {ranked, topIndices, taxaOrder, functionOrder} = heatmapLinks(indices);
  // Full lineages identify rows; short labels must not merge taxa with identical names.
  const axisLabels = Object.fromEntries(ranked.map(index => current.pairs[index])
    .map(pair => [pair.taxon, pairLabel(pair)]));
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
      top_link: topIndices.has(index)
    });
  }
  for (const taxon of taxaOrder) {
    for (const func of functionOrder) {
      if (!seen.get(taxon).has(func)) {
        cells.push({taxon, function: func, present: false, top_link: false});
      }
    }
  }
  // Tooltip content and the selected load are data supplied to the template.
  for (const cell of cells) {
    cell.load = cell.present ? cell[field] : null;
    cell.tooltip = {
      "Taxon": cell.taxon,
      "Function": cell.function,
      [`${label} load`]: cell.present ? cell.load.toExponential(2) : "No estimated load",
      "Top N link": cell.top_link ? "Yes" : "No"
    };
  }
  const spec = plotSpec("heatmap", cells);
  // Bind the same category order and labels to every layer sharing the axes.
  for (const encoding of [spec.encoding, ...spec.layer.map(layer => layer.encoding).filter(encoding => encoding.x && encoding.y)]) {
    encoding.x.sort = functionOrder;
    encoding.y.sort = taxaOrder;
    encoding.y.axis.labelLimit = Math.floor(spec.width * 0.3);
    encoding.y.axis.labelExpr = `${JSON.stringify(axisLabels)}[datum.label]`;
    encoding.x.axis.labels = functionOrder.length <= 30;
  }
  spec.layer[0].encoding.color.title = `${label} load`;
  return spec;
}

/** Clone templates so changing filters never mutates the loaded definitions. */
function plotSpec(name, values) {
  const spec = JSON.parse(JSON.stringify(plotTemplates[name]));
  spec.data.values = values;
  spec.width = plotWidth(name);
  return spec;
}

function sampleGroupDomain(samples) {
  // Include categories omitted from the top 100 bars to keep colors consistent.
  return [...new Set(samples.map(sample => sample.group))].sort();
}

function boxplotSpec(samples) {
  const spec = plotSpec("boxplot", samples);
  const domain = sampleGroupDomain(samples);
  spec.encoding.x.title = groupSelect.value || "Samples";
  spec.encoding.x.sort = domain;
  spec.encoding.color.scale.domain = domain;
  return spec;
}

function sampleSpec(samples) {
  const top = [...samples].sort((a, b) => b.load - a.load).slice(0, 100);
  const spec = plotSpec("samples", top);
  spec.encoding.x.axis.labels = top.length <= 30;
  spec.encoding.color.scale.domain = sampleGroupDomain(samples);
  spec.encoding.tooltip.find(field => field.field === "group").title = groupSelect.value || "Group";
  return spec;
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
    `Each square is one pair at the selected taxonomy level; color shows its ${method} load across all samples, including zeros. The top links are ranked before taxon/function dropdown selections; those selections narrow the axes. All intersections between the remaining axes are shown. Dots in the upper-right corners mark the top links, white cells have a zero summary, and gray cells have no estimated load in the TFA table. Taxa are summed within each sample before this summary. Click a taxon or function axis label to select it, or click a square to inspect that pair.`;
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
      if (item && item.mark && item.mark.role === "axis-label") {
        const scale = item.mark.group.datum.scale;
        const select = scale === "x" ? functionSelect : scale === "y" ? taxonSelect : null;
        // Use the underlying category, not the shortened or truncated label.
        const value = String(item.datum.value);
        if (select && [...select.options].some(option => option.value === value)) {
          select.value = value;
          select.dispatchEvent(new Event("change"));
        }
        return;
      }
      const pair = item && item.datum;
      if (pair && pair.present && pair.taxon && pair.function) {
        taxonSearch.value = "";
        functionSearch.value = "";
        functionSelect.value = pair.function;
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
  refreshSelectors();
  showSelectedTaxonomy();
  showSelectedFunction();
  // Serialize asynchronous redraws so an older render cannot overwrite a newer one.
  renderQueue = renderQueue.then(draw).catch(error => {
    status.textContent = `Could not draw the plots: ${error.message}`;
  });
}

clearSelections.addEventListener("click", resetSelections);
// Rebuild both dropdowns and plots whenever a filter or display limit changes.
showUnclassified.addEventListener("change", () => {
  update();
});
levelSelect.addEventListener("change", () => {
  changeLevel();
  update();
});
taxonSearch.addEventListener("input", () => {
  update();
});
functionSearch.addEventListener("input", () => {
  update();
});
functionSelect.addEventListener("change", () => {
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
