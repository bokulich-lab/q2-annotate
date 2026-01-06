$(document).ready(function () {
    removeBS3refs();
    // adjustTagsToBS3();

    // Data injected from Python
    const taxa = window.taxa || [];
    const samples = window.samples || [];
    const vegaAbundanceHistogramSpec = window.vegaAbundanceHistogramSpec || {};

    const taxonDropdown = document.getElementById('taxonDropdown');
    const sampleDropdown = document.getElementById('sampleDropdown');
    const histogramContainer = document.getElementById('vega-histogram');
    const statisticsDisplay = document.getElementById('statistics-display');
    let vegaView = null;
    let abundanceData = null; // Cache the loaded data

    // Populate dropdowns with "All" option
    const allTaxonOpt = document.createElement('option');
    allTaxonOpt.value = '';
    allTaxonOpt.textContent = 'All';
    taxonDropdown.appendChild(allTaxonOpt);
    
    if (taxa && Array.isArray(taxa) && taxa.length > 0) {
        taxa.forEach(taxon => {
            const opt = document.createElement('option');
            opt.value = taxon;
            opt.textContent = taxon;
            taxonDropdown.appendChild(opt);
        });
    } else {
        console.warn('No taxa data available for dropdown');
    }

    const allSampleOpt = document.createElement('option');
    allSampleOpt.value = '';
    allSampleOpt.textContent = 'All';
    sampleDropdown.appendChild(allSampleOpt);
    
    if (samples && Array.isArray(samples) && samples.length > 0) {
        samples.forEach(sample => {
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

    function calculateStatistics(filteredData) {
        if (!filteredData || filteredData.length === 0) {
            return { count: 0, mean: 0, median: 0, std: 0 };
        }
        
        const abundances = filteredData.map(d => d.abundance).filter(v => v != null && !isNaN(v));
        
        if (abundances.length === 0) {
            return { count: 0, mean: 0, median: 0, std: 0 };
        }
        
        const count = abundances.length;
        const mean = d3.mean(abundances);
        const median = d3.median(abundances);
        const std = d3.deviation(abundances) || 0;
        
        return { count, mean, median, std };
    }

    function updateStatistics(filteredData) {
        const stats = calculateStatistics(filteredData);
        if (stats.count === 0) {
            statisticsDisplay.innerHTML = '<p class="text-muted">No data available for selected filters.</p>';
            return;
        }
        
        statisticsDisplay.innerHTML = `
            <p><strong>Count:</strong> <span id="stat-count">${stats.count.toLocaleString()}</span></p>
            <p><strong>Mean:</strong> <span id="stat-mean">${stats.mean.toFixed(4)}</span></p>
            <p><strong>Median:</strong> <span id="stat-median">${stats.median.toFixed(4)}</span></p>
            <p><strong>Std Dev:</strong> <span id="stat-std">${stats.std.toFixed(4)}</span></p>
        `;
    }

    function updateVisualization() {
        const selectedTaxon = taxonDropdown.value;
        const selectedSample = sampleDropdown.value;

        if (!vegaView) {
            return;
        }

        // Update Vega signals (use empty string if not selected to show all)
        vegaView.signal('taxon_param', selectedTaxon || '')
            .signal('sample_param', selectedSample || '')
            .runAsync();

        // Calculate statistics from filtered data
        // Use cached data loaded directly from Arrow file to avoid Vega dataset access errors
        if (!abundanceData) {
            statisticsDisplay.innerHTML = '<p class="text-muted">Data not yet loaded.</p>';
            return;
        }
        
        // Flatten the abundances arrays and filter based on selections
        let filteredData = [];
        for (const row of abundanceData) {
            // Check if this row matches the filters
            const matchesTaxon = !selectedTaxon || row.taxon === selectedTaxon;
            const matchesSample = !selectedSample || row.sample === selectedSample;
            
            if (matchesTaxon && matchesSample && row.abundances) {
                // Handle both array and typed array formats from Arrow
                let abundances = row.abundances;
                if (!Array.isArray(abundances)) {
                    // Convert typed array to regular array if needed
                    abundances = Array.from(abundances);
                }
                
                // Flatten the abundances array
                for (const abundance of abundances) {
                    const val = typeof abundance === 'number' ? abundance : parseFloat(abundance);
                    if (!isNaN(val) && val > 0) { // Only include valid non-zero abundances
                        filteredData.push({
                            taxon: row.taxon,
                            sample: row.sample,
                            abundance: val
                        });
                    }
                }
            }
        }
        
        if (filteredData.length === 0) {
            statisticsDisplay.innerHTML = '<p class="text-muted">No data available for selected filters.</p>';
            return;
        }
        
        updateStatistics(filteredData);
    }

    const debouncedUpdateVisualization = debounce(updateVisualization, 150);

    // Event listeners for dropdowns
    taxonDropdown.addEventListener('change', debouncedUpdateVisualization);
    sampleDropdown.addEventListener('change', debouncedUpdateVisualization);

    // Debug: Log data availability
    console.log('Taxa:', taxa);
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
            
            // Initialize Vega-Lite visualization after data is loaded
            if (Object.keys(vegaAbundanceHistogramSpec).length > 0) {
                // Validate spec before embedding
                try {
                    vegaEmbed('#vega-histogram', vegaAbundanceHistogramSpec, {
                        actions: true,
                        renderer: 'canvas'
                    }).then(res => {
                        vegaView = res.view;
                        document.getElementById('spinner-histogram')?.remove();
                        
                        // Initial update with loaded data
                        updateVisualization();
                    }).catch(error => {
                        console.error('Error embedding Vega view:', error);
                        console.error('Vega spec:', JSON.stringify(vegaAbundanceHistogramSpec, null, 2));
                        document.getElementById('spinner-histogram')?.remove();
                        histogramContainer.innerHTML = '<p class="text-danger">Error loading visualization: ' + (error.message || error) + '</p>';
                    });
                } catch (error) {
                    console.error('Error preparing Vega spec:', error);
                    document.getElementById('spinner-histogram')?.remove();
                    histogramContainer.innerHTML = '<p class="text-danger">Error: Invalid visualization specification.</p>';
                }
            } else {
                console.error('Vega spec is empty or invalid');
                document.getElementById('spinner-histogram')?.remove();
                histogramContainer.innerHTML = '<p class="text-danger">Error: Visualization specification is empty or invalid.</p>';
            }
        })
        .catch(error => {
            console.error('Error loading abundance data:', error);
            document.getElementById('spinner-histogram')?.remove();
            histogramContainer.innerHTML = '<p class="text-danger">Error loading data: ' + (error.message || error) + '</p>';
        });
});

