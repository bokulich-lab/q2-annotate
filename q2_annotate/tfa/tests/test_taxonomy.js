"use strict";
const assert = require("node:assert/strict");
const {test} = require("node:test");
const {collapseTaxonomy, latestTaxonLabel, isUnclassifiedTaxon} = require("../../assets/tfa_explore/taxonomy.js");
const data = require("./data/taxonomy-collapse.json");
// Normalize structures returned from the isolated browser context.
const plain = value => JSON.parse(JSON.stringify(value));
// Compare original pair summaries independently of heatmap display annotations.
const originalPairs = values => plain(values.filter(pair => pair.present).map(
  ({present, top_link, load, tooltip, ...pair}) => pair
));
const expectedCollapse = require("./data/expected-taxonomy-collapse.json");

test("original taxa preserve the input data", () => {
  assert.equal(collapseTaxonomy(data, 0), data);
});

test("collapse sums within samples and preserves parents, functions, and missing ranks", () => {
  assert.deepEqual(collapseTaxonomy(data, 3), expectedCollapse["3"]);
});

test("total load is conserved at each depth and zeros are retained in summaries", () => {
  for (let level = 1; level <= 4; level++) {
    const result = collapseTaxonomy(data, level);
    assert.equal(result.pairs.reduce((total, pair) => total + pair.total, 0), 22);
    assert.equal(result.loads.reduce((total, load) => total + load[2], 0), 22);
  }
  assert.deepEqual(collapseTaxonomy(data, 1), expectedCollapse["1"]);
});

test("changing the level refreshes the dropdown and every plot, including heatmap clicks", async () => {
  const fs = require("node:fs");
  const path = require("node:path");
  const vm = require("node:vm");
  const elements = new Map();
  const specs = new Map();
  const clicks = new Map();
  const element = id => {
    if (!elements.has(id)) elements.set(id, {
      value: ({level: "0", metric: "sum", limit: "20"})[id] || "",
      textContent: "",
      checked: true, clientWidth: 1000, parentElement: {clientWidth: 1000}, options: [], listeners: {},
      add(option) { this.options.push(option); },
      appendChild(option) { this.options.push(option); },
      replaceChildren() { this.options = []; this.value = ""; },
      addEventListener(event, fn) { this.listeners[event] = fn; }
    });
    return elements.get(id);
  };
  const context = vm.createContext({
    fetch: async url => {
      if (url === "plots.json") return {
        ok: true, json: async () => require("../../assets/tfa_explore/plots.json")
      };
      assert.equal(url, "data.json");
      return {ok:true, json:async () => data};
    },
    document: {getElementById: element, createElement: () => ({})},
    window: {addEventListener() {}},
    Option: function(text, value) { this.textContent = text; this.value = value; },
    collapseTaxonomy, latestTaxonLabel, isUnclassifiedTaxon,
    vegaEmbed: async (id, spec) => {
      specs.set(id, spec);
      return {view: {finalize() {}, addEventListener(event, fn) { clicks.set(id, fn); }}};
    }
  });
  vm.runInContext(fs.readFileSync(path.join(__dirname, "../../assets/tfa_explore/explore.js"), "utf8"), context);
  await vm.runInContext("initialization", context);
  const viewState = () => plain({
    controls: Object.fromEntries(["level", "taxon-search", "taxon", "function", "group", "metric", "limit"]
      .map(id => [id, element(id).value])),
    showUnclassified: element("show-unclassified").checked,
    resetDisabled: element("clear-selections").disabled,
    taxa: element("taxon").options,
    text: Object.fromEntries(["selected-taxonomy", "status", "heatmap-caption", "group-caption"]
      .map(id => [id, element(id).textContent])),
    plots: Object.fromEntries(specs)
  });
  const initialView = viewState();
  assert.deepEqual(plain({
    levels:element("level").options.map(option => [option.textContent, option.value]),
    selected:element("level").value,
    showUnclassified:element("show-unclassified").checked,
    labels:element("taxon").options.map(option => option.textContent),
    heatmap:specs.get("#heatmap-plot").data.values.filter(pair => pair.present).map(pair => ({taxon:pair.taxon_label, function:pair.function, total:pair.total}))
  }), {
    levels:[["Level 1", "1"], ["Level 2", "2"], ["Level 3", "3"], ["Level 4", "4"]],
    selected:"4", showUnclassified:false,
    labels:["All taxa", "s__one", "s__two"],
    heatmap:[
      {taxon:"s__two", function:"geneA", total:5},
      {taxon:"s__one", function:"geneA", total:3},
      {taxon:"s__one", function:"geneB", total:2}
    ]
  });
  // Top links select axes, while below-cutoff values remain visible and clickable.
  const fullHeatmap = specs.get("#heatmap-plot");
  assert.equal(fullHeatmap.data.values.length, 4);
  const absent = fullHeatmap.data.values.find(pair => !pair.present);
  assert.equal(absent.tooltip["Sum load"], "No estimated load");
  clicks.get("#heatmap-plot")({}, {datum: absent});
  assert.equal(element("taxon").value, "");
  assert.equal(element("function").value, "");
  element("limit").value = "2";
  element("limit").listeners.change();
  await vm.runInContext("renderQueue", context);
  // The two highest links both use geneA, so geneB is outside the chosen axes.
  assert.equal(specs.get("#heatmap-plot").data.values.length, 2);
  element("limit").value = "20";
  element("limit").listeners.change();
  await vm.runInContext("renderQueue", context);

  element("level").value = "3";
  element("level").listeners.change();
  await vm.runInContext("renderQueue", context);
  const expectedTaxa = ["", ...[
    "d__Bacteria; p__A; Unclassified (level 3)", "d__Bacteria; p__A; g__Shared",
    "d__Bacteria; p__B; Unclassified (level 3)", "d__Bacteria; p__B; g__Shared"
  ]];
  element("show-unclassified").checked = true;
  element("show-unclassified").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), expectedTaxa);
  element("show-unclassified").checked = false;
  element("show-unclassified").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), [
    "", "d__Bacteria; p__A; g__Shared", "d__Bacteria; p__B; g__Shared"
  ]);
  assert.deepEqual(originalPairs(specs.get("#heatmap-plot").data.values), [
    expectedCollapse["3"].pairs[0], expectedCollapse["3"].pairs[2], expectedCollapse["3"].pairs[1]
  ]);
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#boxplot").data.values)), [
    {sample:"S1", group:"All samples", load:10},
    {sample:"S2", group:"All samples", load:5},
    {sample:"S3", group:"All samples", load:2}
  ]);
  element("show-unclassified").checked = true;
  element("show-unclassified").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), expectedTaxa);
  // Function filtering clears a taxon that no longer belongs in the list.
  element("taxon").value = "d__Bacteria; p__B; g__Shared";
  element("function").value = "geneB";
  element("function").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual({
    selected:element("taxon").value,
    options:element("taxon").options.map(option => ({value:option.value, label:option.textContent}))
  }, {
    selected:"", options:[
      {value:"", label:"All taxa"},
      {value:"d__Bacteria; p__A; g__Shared", label:"g__Shared"}
    ]
  });
  // Search filters all plots by full lineage, irrespective of case and whitespace.
  element("function").value = "";
  element("function").listeners.change();
  element("taxon-search").value = "  P__b  ";
  element("taxon-search").listeners.input();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), [
    "", "d__Bacteria; p__B; Unclassified (level 3)", "d__Bacteria; p__B; g__Shared"
  ]);
  assert.equal(element("taxon").options[0].textContent, "All matching taxa");
  assert.deepEqual(originalPairs(specs.get("#heatmap-plot").data.values), [
    expectedCollapse["3"].pairs[2], expectedCollapse["3"].pairs[4]
  ]);
  assert.deepEqual(plain(specs.get("#boxplot").data.values).map(sample => sample.load), [7, 4, 0]);
  assert.deepEqual(plain(specs.get("#samples").data.values).map(sample => sample.load), [7, 4, 0]);
  // Rank within search results: the globally highest link is outside p__B.
  element("limit").value = "1";
  element("limit").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(originalPairs(specs.get("#heatmap-plot").data.values), [expectedCollapse["3"].pairs[2]]);
  element("limit").value = "20";
  element("limit").listeners.change();
  element("function").value = "geneB";
  element("function").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), [""]);
  assert.equal(element("heatmap-plot").textContent, "No values to plot for this selection.");
  element("function").value = "";
  element("function").listeners.change();
  element("taxon-search").value = "no matching taxon";
  element("taxon-search").listeners.input();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(element("taxon").options.map(option => option.value), [""]);
  for (const id of ["heatmap-plot", "boxplot", "samples"]) {
    assert.equal(element(id).textContent, "No values to plot for this selection.");
  }
  // Clearing search restores all taxa and permits drilling into a heatmap cell.
  element("taxon-search").value = "";
  element("taxon-search").listeners.input();
  await vm.runInContext("renderQueue", context);
  assert.equal(specs.get("#heatmap-plot").data.values.filter(pair => pair.present).length, 5);
  const shared = specs.get("#heatmap-plot").data.values.find(pair => pair.taxon === "d__Bacteria; p__A; g__Shared" && pair.function === "geneA");
  clicks.get("#heatmap-plot")({}, {datum: shared});
  await vm.runInContext("renderQueue", context);
  assert.deepEqual({
    taxon:element("taxon").value, search:element("taxon-search").value,
    lineage:element("selected-taxonomy").textContent
  }, {taxon:shared.taxon, search:"", lineage:`Selected lineage: ${shared.taxon}`});
  assert.equal(specs.get("#heatmap-plot").encoding.y.field, "taxon");
  assert.equal(specs.get("#heatmap-plot").encoding.y.axis.labelLimit, 300);
  const expectedSamples = [
    {sample:"S1", group:"All samples", load:3},
    {sample:"S2", group:"All samples", load:5},
    {sample:"S3", group:"All samples", load:0}
  ];
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#boxplot").data.values)), expectedSamples);
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#samples").data.values)), [
    expectedSamples[1], expectedSamples[0], expectedSamples[2]
  ]);
  element("metric").value = "median";
  element("metric").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#heatmap-plot").layer[0].encoding.color)), {
    condition:[{test:"!datum.present", value:"#f3f3f3"}, {test:"datum.load === 0", value:"#ffffff"}],
    field:"load", type:"quantitative", title:"Median load", scale:{scheme:"blues", zero:true}
  });
  assert.deepEqual(originalPairs(specs.get("#heatmap-plot").data.values), [expectedCollapse["3"].pairs[0]]);
  element("function").value = "";
  element("function").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(originalPairs(specs.get("#heatmap-plot").data.values), expectedCollapse["3"].pairs.slice(0,2));
  element("level").value = "1";
  element("level").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual({selected:element("taxon").value, options:element("taxon").options.map(option => option.value)}, {
    selected:"", options:["", "d__Bacteria"]
  });
  element("limit").value = "1000";
  element("limit").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(originalPairs(specs.get("#heatmap-plot").data.values), expectedCollapse["1"].pairs);
  assert.deepEqual(JSON.parse(JSON.stringify(specs.get("#boxplot").data.values)), [
    {sample:"S1", group:"All samples", load:11},
    {sample:"S2", group:"All samples", load:9},
    {sample:"S3", group:"All samples", load:2}
  ]);
  element("group").value = "treatment";
  element("group").listeners.change();
  await vm.runInContext("renderQueue", context);
  const boxColor = specs.get("#boxplot").encoding.color;
  const barColor = specs.get("#samples").encoding.color;
  const expectedColor = {
    field:"group", type:"nominal", legend:null,
    scale:{domain:["A", "B"], scheme:"tableau10"}
  };
  assert.deepEqual(JSON.parse(JSON.stringify({box:boxColor, bars:barColor})), {
    box:expectedColor, bars:expectedColor
  });
  // A low-load category omitted from the top 100 must not shift other colors.
  const limited = JSON.parse(JSON.stringify(vm.runInContext(`(() => {
    const samples = Array.from({length:101}, (_, index) => ({
      sample: String(index), group: index === 100 ? "Absent from bars" : "Visible",
      load: 101 - index
    }));
    return {box: boxplotSpec(samples), bars: sampleSpec(samples)};
  })()`, context)));
  const expectedLimitedColor = {
    field:"group", type:"nominal", legend:null,
    scale:{domain:["Absent from bars", "Visible"], scheme:"tableau10"}
  };
  assert.deepEqual({box:limited.box.encoding.color, bars:limited.bars.encoding.color}, {
    box:expectedLimitedColor, bars:expectedLimitedColor
  });
  assert.deepEqual(limited.bars.data.values, Array.from({length:100}, (_, index) => ({
    sample:String(index), group:"Visible", load:101-index
  })));

  // Reset filters and display settings together, including a selected taxon.
  element("function").value = "geneA";
  element("function").listeners.change();
  element("taxon-search").value = "bacteria";
  element("taxon-search").listeners.input();
  element("taxon").value = "d__Bacteria";
  element("taxon").listeners.change();
  element("show-unclassified").checked = true;
  element("show-unclassified").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.notDeepEqual(viewState(), initialView);
  element("clear-selections").listeners.click();
  await vm.runInContext("renderQueue", context);
  assert.deepEqual(viewState(), initialView);

  // A diagonal top two must retain the third (off-diagonal) positive pair.
  const matrix = plain(vm.runInContext(`(() => {
    const original = current, originalLoads = loadsByPair;
    const originalMetric = metricSelect.value, originalLimit = limitSelect.value;
    try {
      current = {pairs: [
        {taxon: "A", taxon_label: "A", function: "X", total: 10, mean: 10/3, median: 0},
        {taxon: "B", taxon_label: "B", function: "Y", total: 9, mean: 3, median: 2},
        {taxon: "A", taxon_label: "A", function: "Y", total: 2, mean: 2/3, median: 0},
        {taxon: "C", taxon_label: "C", function: "Z", total: 0.000123456, mean: 0.000123456/3, median: 0}
      ]};
      loadsByPair = [[[0, 0, 10]], [[1, 0, 7], [1, 1, 2]], [[2, 0, 2]], [[3, 0, 0.000123456]]];
      limitSelect.value = "2";
      metricSelect.value = "sum";
      const sum = heatmapSpec([0, 1, 2, 3]);
      const small = heatmapSpec([3]);
      metricSelect.value = "median";
      const median = heatmapSpec([0, 1, 2, 3]);
      const filtered = heatmapSpec([2]);
      limitSelect.value = "1000";
      const all = heatmapSpec([0, 1, 2, 3]);
      return {sum, small, median, filtered, all};
    } finally {
      current = original; loadsByPair = originalLoads;
      metricSelect.value = originalMetric; limitSelect.value = originalLimit;
    }
  })()`, context));
  const cell = (spec, taxon, func) => spec.data.values.find(
    pair => pair.taxon === taxon && pair.function === func
  );
  assert.equal(matrix.sum.data.values.length, 4);
  assert.equal(cell(matrix.sum, "A", "X").tooltip["Sum load"], "1.00e+1");
  assert.equal(cell(matrix.small, "C", "Z").tooltip["Sum load"], "1.23e-4");
  assert.equal(cell(matrix.small, "C", "Z").load, 0.000123456);
  assert.equal(cell(matrix.sum, "A", "Y").total, 2);
  assert.equal(cell(matrix.sum, "A", "Y").top_link, false);
  assert.equal(cell(matrix.sum, "B", "X").present, false);
  assert.equal(cell(matrix.sum, "B", "X").total, undefined);
  assert.equal(cell(matrix.sum, "C", "Z"), undefined);
  assert.deepEqual(matrix.sum.data.values.filter(pair => pair.top_link).map(
    pair => [pair.taxon, pair.function]), [["A", "X"], ["B", "Y"]]);
  assert.equal(cell(matrix.median, "A", "X").load, 0);
  assert.equal(cell(matrix.median, "A", "X").tooltip["Median load"], "0.00e+0");
  assert.equal(cell(matrix.median, "A", "X").total, 10);
  assert.equal(matrix.filtered.data.values.length, 1);
  assert.equal(cell(matrix.filtered, "A", "Y").present, true);
  assert.equal(matrix.all.data.values.length, 9);
  assert.equal(matrix.all.data.values.filter(pair => pair.top_link).length, 4);

  // Count displayed functions, rather than pairs, when hiding crowded labels.
  const heatmapAxes = plain(vm.runInContext(`(() => {
    const original = current;
    const originalLimit = limitSelect.value;
    try {
      return [[30, "1000"], [31, "1000"], [31, "20"]].map(([count, limit]) => {
        current = {pairs: Array.from({length:count * 2}, (_, index) => ({
          ...original.pairs[0], taxon: "taxon" + index % 2,
          function: "gene" + Math.floor(index / 2), total: count - Math.floor(index / 2)
        }))};
        limitSelect.value = limit;
        return heatmapSpec(current.pairs.map((_, index) => index)).encoding.x;
      });
    } finally {
      current = original;
      limitSelect.value = originalLimit;
    }
  })()`, context));
  assert.deepEqual(heatmapAxes, [30, 31, 10].map(count => ({
    field:"function", type:"nominal", sort:Array.from({length:count}, (_, index) => "gene" + index),
    title:"Function", scale:{paddingInner:0, paddingOuter:0}, axis:{labels:count <= 30}
  })));
});


test("deep levels use the last classified rank as the heatmap label", () => {
  const result = collapseTaxonomy(data, 10);
  assert.deepEqual(result.pairs.map(pair => ({id:pair.taxon_id, gene:pair.function, label:pair.taxon_label})), [
    {id:"T1", gene:"geneA", label:"s__one"}, {id:"T1", gene:"geneB", label:"s__one"},
    {id:"T2", gene:"geneA", label:"s__two"}, {id:"T3", gene:"geneA", label:"g__Shared"},
    {id:"T4", gene:"geneA", label:"p__A"}, {id:"T5", gene:"geneA", label:"p__B"}
  ]);
  assert.deepEqual([
    latestTaxonLabel(["d__Bacteria", "g__", "s__unclassified"], "ID"),
    latestTaxonLabel(["", "g__"], "ID"),
    latestTaxonLabel(["d__Unclassified", "", ""], "ID"),
    latestTaxonLabel(["Unclassified"], "ID"),
    latestTaxonLabel(["ID"], "ID")
  ], ["d__Bacteria", "Unclassified", "d__Unclassified", "Unclassified", "Unclassified"]);
});


test("unclassified taxonomy keeps its label when collapsed to a deeper level", () => {
  const unclassified = require("./data/unclassified-taxonomy.json");
  assert.deepEqual(collapseTaxonomy(unclassified, 3), {
    samples:["S1", "S2"], groups:{},
    pairs:[{
      taxon:"d__Unclassified; Unclassified (level 2); Unclassified (level 3)",
      unclassified:true, taxon_label:"d__Unclassified", taxon_short:"Unclassified (level 3)",
      taxon_id:"T013", function:"geneA", total:3, mean:1.5, median:1.5
    }],
    loads:[[0, 0, 3]]
  });
});


test("unclassified includes ID-only taxa and taxa missing the selected rank", () => {
  assert.deepEqual([
    isUnclassifiedTaxon({lineage:["T013"], taxon_id:"T013"}, 0),
    isUnclassifiedTaxon(data.pairs[4], 0), isUnclassifiedTaxon(data.pairs[4], 3),
    isUnclassifiedTaxon(data.pairs[5], 3), isUnclassifiedTaxon(data.pairs[0], 3),
    isUnclassifiedTaxon({lineage:["d__Bacteria", "g__unclassified"]}, 0)
  ], [true, false, true, true, false, true]);
});

test("JSON loading failures show a message and keep controls disabled", async () => {
  const fs = require("node:fs");
  const path = require("node:path");
  const vm = require("node:vm");
  const source = fs.readFileSync(path.join(__dirname, "../../assets/tfa_explore/explore.js"), "utf8");
  for (const [fetch, message] of [
    [async () => ({ok:false, status:404}), "data.json: HTTP 404"],
    [async url => url === "plots.json"
      ? {ok:false, status:404} : {ok:true, json:async () => data}, "plots.json: HTTP 404"],
    [async url => ({ok:true, json:async () => {
      if (url === "plots.json") throw new Error("Invalid plot JSON");
      return data;
    }}), "Invalid plot JSON"],
    [async () => { throw new Error("Network unavailable"); }, "Network unavailable"],
    [async () => ({ok:true, json:async () => { throw new Error("Invalid JSON"); }}), "Invalid JSON"]
  ]) {
    const elements = new Map();
    const getElementById = id => {
      if (!elements.has(id)) elements.set(id, {value:"", addEventListener() {}});
      return elements.get(id);
    };
    const context = vm.createContext({fetch, document:{getElementById}, window:{addEventListener() {}}});
    vm.runInContext(source, context);
    await vm.runInContext("initialization", context);
    assert.equal(getElementById("status").textContent, `Could not load TFA data: ${message}`);
    assert.deepEqual(Object.fromEntries([
      "show-unclassified", "level", "taxon-search", "taxon", "function", "group", "metric", "limit", "clear-selections"
    ].map(id => [id, getElementById(id).disabled])), {
      "show-unclassified":true, level:true, "taxon-search":true, taxon:true,
      function:true, group:true, metric:true, limit:true, "clear-selections":true
    });
  }
});
