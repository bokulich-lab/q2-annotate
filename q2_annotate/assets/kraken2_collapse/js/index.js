$(document).ready(function () {
    removeBS3refs();

    // Data injected from Python
    const samples = window.samples || [];
    const vegaAbundanceHistogramSpec = window.vegaAbundanceHistogramSpec || {};

    const sampleDropdown = document.getElementById('sampleDropdown');
    const taxonDropdown = document.getElementById('taxonDropdown');
    const histogramGrid = document.getElementById('histogram-grid');
    const selectionStats = document.getElementById('selection-stats');
    let vegaViews = {}; // Store views for each taxon
    let abundanceData = null; // Cache the loaded data
    let taxaList = []; // List of unique taxa
    let taxaMap = {}; // Map from taxon short name to full taxon name

    // Populate sample dropdown with "All" option
    const allSampleOpt = document.createElement('option');
    allSampleOpt.value = '';
    allSampleOpt.textContent = 'All';
    sampleDropdown.appendChild(allSampleOpt);
    
    if (samples && Array.isArray(samples) && samples.length > 0) {
        // Sort samples alphabetically
        const sortedSamples = [...samples].sort();
        sortedSamples.forEach(sample => {
            const opt = document.createElement('option');
            opt.value = sample;
            opt.textContent = sample;
            sampleDropdown.appendChild(opt);
        });
    } else {
        console.warn('No samples data available for dropdown');
    }

    function debounce(func, delay) {
        let timeout;
        return function(...args) {
            const context = this;
            clearTimeout(timeout);
            timeout = setTimeout(() => func.apply(context, args), delay);
        };
    }

    function calculateStatisticsForSelection(taxon, sample) {
        // Get filtered data for the selected taxon and sample
        if (!abundanceData || abundanceData.length === 0) {
            return null;
        }
        
        let filteredData = [];
        for (const row of abundanceData) {
            const matchesTaxon = !taxon || row.taxon === taxon;
            const matchesSample = !sample || row.sample === sample;
            
            if (matchesTaxon && matchesSample && row.abundances) {
                let abundances = row.abundances;
                if (!Array.isArray(abundances)) {
                    abundances = Array.from(abundances);
                }
                
                for (const abundance of abundances) {
                    const val = typeof abundance === 'number' ? abundance : parseFloat(abundance);
                    if (!isNaN(val) && val > 0) {
                        filteredData.push(val);
                    }
                }
            }
        }
        
        if (filteredData.length === 0) {
            return null;
        }
        
        const count = filteredData.length;
        const mean = d3.mean(filteredData);
        const median = d3.median(filteredData);
        const std = d3.deviation(filteredData) || 0;
        
        return { count, mean, median, std };
    }

    function updateSelectionStats() {
        const selectedTaxonLabel = taxonDropdown.value;
        const selectedSample = sampleDropdown.value || '';
        const selectedTaxon = selectedTaxonLabel ? (taxaMap[selectedTaxonLabel] || selectedTaxonLabel) : '';
        
        const stats = calculateStatisticsForSelection(selectedTaxon, selectedSample);
        
        if (!stats) {
            selectionStats.innerHTML = `
                <div class="stats-title">Abundance metrics</div>
                <small class="text-muted">Select taxon and sample to view statistics</small>
            `;
            return;
        }
        
        selectionStats.innerHTML = `
            <div class="stats-title">Abundance metrics</div>
            <div class="stats-content">
                <div class="stat-item"><strong>Contig count:</strong> <span>${stats.count.toLocaleString()}</span></div>
                <div class="stat-item"><strong>Mean:</strong> <span>${stats.mean.toFixed(1)}</span></div>
                <div class="stat-item"><strong>Median:</strong> <span>${stats.median.toFixed(1)}</span></div>
                <div class="stat-item"><strong>Std Dev:</strong> <span>${stats.std.toFixed(1)}</span></div>
            </div>
        `;
    }

    function getTaxonShort(taxon) {
        // Extract last taxonomy level (species) without prefix for plot titles
        const parts = taxon.split(';');
        let shortName = parts.length > 0 ? parts[parts.length - 1] : taxon;
        
        // Strip prefix by splitting on "__" and taking the right side
        if (shortName.includes('__')) {
            shortName = shortName.split('__').slice(-1)[0];
        }
        
        return shortName;
    }

    function getTaxonDropdownLabel(taxon) {
        // Extract last two taxonomy levels with their level descriptors (e.g., "g__Genus;s__Species")
        const parts = taxon.split(';');
        if (parts.length >= 2) {
            // Return last two levels joined by semicolon
            return parts.slice(-2).join(';');
        }
        // If only one level, return it as-is
        return parts.length > 0 ? parts[parts.length - 1] : taxon;
    }

    function hasDataForTaxon(taxon, sampleFilter = '') {
        // Check if there's any data for this taxon (optionally filtered by sample)
        if (!abundanceData || abundanceData.length === 0) {
            return false;
        }
        
        for (const row of abundanceData) {
            if (row.taxon === taxon) {
                // If no sample filter, any data for this taxon counts
                if (!sampleFilter) {
                    if (row.abundances && (Array.isArray(row.abundances) ? row.abundances.length > 0 : Array.from(row.abundances).length > 0)) {
                        return true;
                    }
                } else {
                    // Check if this row matches the sample filter and has data
                    if (row.sample === sampleFilter && row.abundances) {
                        const abundances = Array.isArray(row.abundances) ? row.abundances : Array.from(row.abundances);
                        if (abundances.length > 0) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    function updateVisualization() {
        const selectedSample = sampleDropdown.value || '';

        // Update all histogram items to show/hide plots or "no data" messages
        const histogramItems = document.querySelectorAll('.histogram-item');
        histogramItems.forEach(item => {
            const taxon = item.dataset.taxon;
            const plotDiv = item.querySelector('.histogram-plot, .histogram-no-data');
            
            if (!taxon || !plotDiv) return;
            
            // Check if there's data for this taxon with the current sample filter
            if (!hasDataForTaxon(taxon, selectedSample)) {
                // Show "no data" message
                if (plotDiv.className !== 'histogram-no-data') {
                    // Remove existing plot if it exists
                    if (vegaViews[taxon]) {
                        delete vegaViews[taxon];
                    }
                    plotDiv.className = 'histogram-no-data';
                    plotDiv.innerHTML = '<p class="text-muted" style="text-align: center; padding: 50px 10px; margin: 0;">No data available</p>';
                }
            } else {
                // There is data - show plot
                if (plotDiv.className === 'histogram-no-data') {
                    // Recreate the plot
                    plotDiv.className = 'histogram-plot';
                    plotDiv.innerHTML = '';
                    
                    const spec = JSON.parse(JSON.stringify(vegaAbundanceHistogramSpec));
                    spec.params[1].value = taxon;
                    
                    vegaEmbed(plotDiv, spec, {
                        actions: false,
                        renderer: 'canvas'
                    }).then(res => {
                        vegaViews[taxon] = res.view;
                        // Set the sample filter
                        res.view.signal('sample_param', selectedSample).runAsync();
                    }).catch(error => {
                        console.error(`Error embedding histogram for taxon ${taxon}:`, error);
                        plotDiv.innerHTML = '<p class="text-danger" style="font-size: 10px;">Error loading plot</p>';
                    });
                } else if (vegaViews[taxon]) {
                    // Update existing plot with new sample filter
                    vegaViews[taxon].signal('sample_param', selectedSample)
                        .runAsync();
                }
            }
        });
    }

    function createHistogramForTaxon(taxon, taxonShort) {
        // Create a container div for this taxon's histogram
        const containerDiv = document.createElement('div');
        containerDiv.className = 'histogram-item';
        containerDiv.dataset.taxon = taxon; // Store taxon for later reference
        
        // Create a title div for the taxon name
        const titleDiv = document.createElement('div');
        titleDiv.className = 'histogram-title';
        titleDiv.textContent = taxonShort;
        containerDiv.appendChild(titleDiv);
        
        // Create a div for the Vega visualization or "no data" message
        const plotDiv = document.createElement('div');
        plotDiv.id = `histogram-${taxon.replace(/[^a-zA-Z0-9]/g, '_')}`;
        plotDiv.className = 'histogram-plot';
        containerDiv.appendChild(plotDiv);
        
        histogramGrid.appendChild(containerDiv);
        
        // Check if there's data for this taxon (with current sample filter)
        const selectedSample = sampleDropdown.value || '';
        if (!hasDataForTaxon(taxon, selectedSample)) {
            plotDiv.className = 'histogram-no-data';
            plotDiv.innerHTML = '<p class="text-muted" style="text-align: center; padding: 50px 10px; margin: 0;">No data available</p>';
            return;
        }
        
        // Create a copy of the spec for this taxon
        const spec = JSON.parse(JSON.stringify(vegaAbundanceHistogramSpec));
        spec.params[1].value = taxon; // Set taxon_param
        
        // Embed the visualization - width: "container" will automatically fit
        vegaEmbed(plotDiv, spec, {
            actions: false,
            renderer: 'canvas'
        }).then(res => {
            vegaViews[taxon] = res.view;
        }).catch(error => {
            console.error(`Error embedding histogram for taxon ${taxon}:`, error);
            plotDiv.innerHTML = '<p class="text-danger" style="font-size: 10px;">Error loading plot</p>';
        });
    }

    function highlightTaxon(selectedTaxon) {
        // Remove highlight from all items
        const histogramItems = document.querySelectorAll('.histogram-item');
        histogramItems.forEach(item => {
            item.classList.remove('highlighted');
        });
        
        // Highlight the selected taxon's container
        if (selectedTaxon) {
            const fullTaxon = taxaMap[selectedTaxon] || selectedTaxon;
            const targetItem = document.querySelector(`.histogram-item[data-taxon="${fullTaxon}"]`);
            if (targetItem) {
                targetItem.classList.add('highlighted');
                // Scroll into view
                targetItem.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
            }
        }
    }

    function updateTaxonHighlight() {
        const selectedTaxon = taxonDropdown.value;
        highlightTaxon(selectedTaxon);
        updateSelectionStats();
    }

    const debouncedUpdateVisualization = debounce(updateVisualization, 150);

    // Event listener for sample dropdown
    sampleDropdown.addEventListener('change', () => {
        debouncedUpdateVisualization();
        updateSelectionStats();
    });
    
    // Event listener for taxon dropdown
    taxonDropdown.addEventListener('change', updateTaxonHighlight);

    // Debug: Log data availability
    console.log('Samples:', samples);
    console.log('Vega spec keys:', Object.keys(vegaAbundanceHistogramSpec));
    
    // Load abundance data directly from Arrow file to avoid Vega dataset access errors
    // This is loaded separately from Vega so we can calculate statistics without triggering Vega errors
    fetch('data/abundance_data.arrow')
        .then(response => {
            if (!response.ok) {
                throw new Error('Failed to load data file');
            }
            return response.arrayBuffer();
        })
        .then(arrayBuffer => {
            // Use vega.read to parse Arrow format
            return vega.read(arrayBuffer, { type: 'arrow', parse: 'auto' });
        })
        .then(data => {
            abundanceData = data;
            console.log('Loaded abundance data:', abundanceData.length, 'rows');
            
            // Extract unique taxa from the data
            const taxaSet = new Set();
            for (const row of abundanceData) {
                if (row.taxon) {
                    taxaSet.add(row.taxon);
                }
            }
            taxaList = Array.from(taxaSet).sort();
            console.log('Found', taxaList.length, 'unique taxa');
            
            // Populate taxon dropdown
            const allTaxonOpt = document.createElement('option');
            allTaxonOpt.value = '';
            allTaxonOpt.textContent = 'None';
            taxonDropdown.appendChild(allTaxonOpt);
            
            // Create map of dropdown labels to full names
            const taxonLabels = [];
            taxaList.forEach(taxon => {
                const taxonShort = getTaxonShort(taxon); // For plot titles
                const taxonLabel = getTaxonDropdownLabel(taxon); // For dropdown (last 2 levels with prefixes)
                taxaMap[taxonLabel] = taxon; // Map dropdown label to full taxon name
                taxonLabels.push(taxonLabel);
            });
            
            // Sort taxon labels alphabetically and populate dropdown
            taxonLabels.sort().forEach(taxonLabel => {
                const opt = document.createElement('option');
                opt.value = taxonLabel;
                opt.textContent = taxonLabel;
                taxonDropdown.appendChild(opt);
            });
            
            // Remove spinner
            document.getElementById('spinner-histogram')?.remove();
            
            // Create a histogram for each taxon
            if (Object.keys(vegaAbundanceHistogramSpec).length > 0 && taxaList.length > 0) {
                taxaList.forEach(taxon => {
                    const taxonShort = getTaxonShort(taxon);
                    createHistogramForTaxon(taxon, taxonShort);
                });
                
                // Initial update with loaded data
                updateVisualization();
                updateSelectionStats();
            } else {
                histogramGrid.innerHTML = '<p class="text-danger">Error: No data or invalid visualization specification.</p>';
            }
        })
        .catch(error => {
            console.error('Error loading abundance data:', error);
            document.getElementById('spinner-histogram')?.remove();
            histogramGrid.innerHTML = '<p class="text-danger">Error loading data: ' + (error.message || error) + '</p>';
        });
});

