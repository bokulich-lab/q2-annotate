"use strict";

// Sum within each sample before calculating summaries across samples.
function collapseTaxonomy(data, level) {
  if (level === 0) return data;
  const groups = new Map();
  const groupByPair = [];
  for (const pair of data.pairs) {
    const lineage = pair.lineage.slice(0, level);
    while (lineage.length < level) lineage.push("");
    const known = lineage.some(label => label && !/^[a-z]__$/i.test(label));
    const key = JSON.stringify([lineage, known ? null : pair.taxon_id]);
    if (!groups.has(key)) {
      const labels = lineage.map((label, index) =>
        label && !/^[a-z]__$/i.test(label) ? label : `Unclassified (level ${index + 1})`);
      const label = labels.join("; ") + (known ? "" : ` [${pair.taxon_id}]`);
      groups.set(key, {taxon: label, taxon_short: labels.at(-1), ids: new Set(), functions: new Map()});
    }
    const group = groups.get(key);
    group.ids.add(pair.taxon_id);
    if (!group.functions.has(pair.function)) {
      group.functions.set(pair.function, new Map());
    }
    groupByPair.push(group.functions.get(pair.function));
  }
  for (const [pair, sample, load] of data.loads) {
    const values = groupByPair[pair];
    values.set(sample, (values.get(sample) || 0) + load);
  }
  const shortCounts = new Map();
  for (const group of groups.values()) {
    shortCounts.set(group.taxon_short, (shortCounts.get(group.taxon_short) || 0) + 1);
  }
  const pairs = [], loads = [];
  for (const group of groups.values()) {
    for (const [func, values] of group.functions) {
      const stored = [...values.values()].filter(value => value !== 0).sort((a, b) => a - b);
      const total = stored.reduce((sum, value) => sum + value, 0);
      const zeros = data.samples.length - stored.length;
      const middle = [(data.samples.length - 1) >> 1, data.samples.length >> 1];
      const median = data.samples.length
        ? middle.reduce((sum, index) => sum + (index < zeros ? 0 : stored[index - zeros]), 0) / 2
        : 0;
      const index = pairs.length;
      pairs.push({
        taxon: group.taxon,
        taxon_short: shortCounts.get(group.taxon_short) > 1 ? group.taxon : group.taxon_short,
        taxon_id: [...group.ids].join(", "),
        function: func,
        total,
        mean: data.samples.length ? total / data.samples.length : 0,
        median
      });
      for (const [sample, load] of values) {
        if (load !== 0) loads.push([index, sample, load]);
      }
    }
  }
  return {...data, pairs, loads};
}

if (typeof module !== "undefined") module.exports = {collapseTaxonomy};
