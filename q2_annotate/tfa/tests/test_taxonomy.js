"use strict";
const assert = require("node:assert/strict");
const {test} = require("node:test");
const {collapseTaxonomy} = require("../../assets/tfa_explore/taxonomy.js");
const data = require("./data/taxonomy-collapse.json");

test("original taxa preserve the input data", () => {
  assert.equal(collapseTaxonomy(data, 0), data);
});

test("collapse sums within samples before computing mean and median", () => {
  const result = collapseTaxonomy(data, 3);
  const index = result.pairs.findIndex(pair => pair.taxon === "d__Bacteria; p__A; g__Shared" && pair.function === "geneA");
  const pair = result.pairs[index];
  assert.equal(pair.total, 8);
  assert.equal(pair.mean, 8 / 3);
  assert.equal(pair.median, 3);
  assert.deepEqual(result.loads.filter(load => load[0] === index), [[index, 0, 3], [index, 1, 5]]);
  assert.equal(result.pairs.filter(pair => pair.function === "geneB").length, 1);
  assert.equal(result.groups, data.groups);
});

test("identical terminal names under different parents remain distinct", () => {
  const pairs = collapseTaxonomy(data, 3).pairs.filter(pair => pair.taxon.endsWith("g__Shared") && pair.function === "geneA");
  assert.equal(pairs.length, 2);
  assert.notEqual(pairs[0].taxon_short, pairs[1].taxon_short);
});

test("missing ranks and empty rank prefixes stay under their known parent", () => {
  const pairs = collapseTaxonomy(data, 3).pairs.filter(pair => pair.taxon.includes("Unclassified"));
  assert.deepEqual(pairs.map(pair => pair.taxon), [
    "d__Bacteria; p__A; Unclassified (level 3)",
    "d__Bacteria; p__B; Unclassified (level 3)"
  ]);
});

test("total load is conserved at each depth and zeros are retained in summaries", () => {
  for (let level = 1; level <= 4; level++) {
    const result = collapseTaxonomy(data, level);
    assert.equal(result.pairs.reduce((total, pair) => total + pair.total, 0), 22);
    assert.equal(result.loads.reduce((total, load) => total + load[2], 0), 22);
  }
  const geneB = collapseTaxonomy(data, 1).pairs.find(pair => pair.function === "geneB");
  assert.equal(geneB.median, 0);
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
      textContent: id === "tfa-data" ? JSON.stringify(data) : "",
      parentElement: {clientWidth: 1000}, options: [], listeners: {},
      add(option) { this.options.push(option); },
      appendChild(option) { this.options.push(option); },
      replaceChildren() { this.options = []; this.value = ""; },
      addEventListener(event, fn) { this.listeners[event] = fn; }
    });
    return elements.get(id);
  };
  const context = vm.createContext({
    document: {getElementById: element, createElement: () => ({})},
    window: {addEventListener() {}},
    Option: function(text, value) { this.textContent = text; this.value = value; },
    collapseTaxonomy,
    vegaEmbed: async (id, spec) => {
      specs.set(id, spec);
      return {view: {finalize() {}, addEventListener(event, fn) { clicks.set(id, fn); }}};
    }
  });
  vm.runInContext(fs.readFileSync(path.join(__dirname, "../../assets/tfa_explore/explore.js"), "utf8"), context);
  await vm.runInContext("renderQueue", context);
  element("level").value = "3";
  element("level").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").options.length, 5); // All taxa plus four lineages.
  const shared = specs.get("#heatmap").data.values.find(pair => pair.taxon === "d__Bacteria; p__A; g__Shared" && pair.function === "geneA");
  clicks.get("#heatmap")({}, {datum: shared});
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").value, shared.taxon);
  assert.deepEqual(Array.from(specs.get("#boxplot").data.values, row => row.load), [3, 5, 0]);
  assert.deepEqual(Array.from(specs.get("#samples").data.values, row => row.load), [5, 3, 0]);
  element("metric").value = "median";
  element("metric").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(specs.get("#heatmap").encoding.color.field, "median");
  assert.equal(specs.get("#heatmap").data.values[0].median, 3);
  element("level").value = "1";
  element("level").listeners.change();
  await vm.runInContext("renderQueue", context);
  assert.equal(element("taxon").value, "");
  assert.equal(element("taxon").options.length, 2);
  assert.deepEqual(Array.from(specs.get("#boxplot").data.values, row => row.load), [11, 9, 0]);
});
