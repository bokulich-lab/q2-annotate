$(document).ready(function () {
    removeBS3refs();

    // ============================================
    // 1. DATA & STATE
    // ============================================

    // Data injected from Python
    const samples = window.samples || [];
    const vegaAbundanceHistogramSpec = window.vegaAbundanceHistogramSpec || {};
    const meanAbundancesData = window.meanAbundances || {};

    // State variables
    let vegaViews = {}; // Store views for each taxon
    let abundanceDataBySample = {}; // Cache the loaded data per sample
    let taxaList = []; // List of unique taxa (sorted by mean abundance)
    let taxaMap = {}; // Map from taxon short name to full taxon name

    // ============================================
    // 2. DOM ELEMENT REFERENCES
    // ============================================

    const sampleDropdown = document.getElementById('sampleDropdown');
    const taxonDropdown = document.getElementById('taxonDropdown');
    const histogramGrid = document.getElementById('histogram-grid');
    const selectionStats = document.getElementById('selection-stats');
    const taxaLimitSlider = document.getElementById('taxaLimitSlider');
    const taxaLimitValue = document.getElementById('taxaLimitValue');

    // ============================================
    // 3. UTILITY FUNCTIONS
    // ============================================

    function debounce(func, delay) {
        let timeout;
        return function(...args) {
            const context = this;
            clearTimeout(timeout);
            timeout = setTimeout(() => func.apply(context, args), delay);
        };
    }

    function getTaxonShort(taxon) {
        // Extract last taxonomy level (species) without prefix for plot titles
        const parts = taxon.split(';');
        const lastPart = parts[parts.length - 1] || taxon;

        // Strip prefix by splitting on "__" and taking the right side
        return lastPart.includes('__') ? lastPart.split('__')[1] : lastPart;
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

    const debouncedUpdateVisualization = debounce(updateVisualization, 150);

    // ============================================
    // 4. DATA LOADING & PROCESSING
    // ============================================

    function loadSampleData(sampleId) {
        // Load JSON data for a sample if not already loaded
        if (abundanceDataBySample[sampleId]) {
            return Promise.resolve();
        }

        const jsonPath = `data/${sampleId}.json`;
        return fetch(jsonPath)
            .then(response => {
                if (!response.ok) {
                    throw new Error(`Failed to load ${jsonPath}`);
                }
                return response.json();
            })
            .then(data => {
                abundanceDataBySample[sampleId] = data;
            });
    }

    function getFilteredAbundanceData(taxon, sample = '') {
        // Helper function to extract filtered abundance values for a given taxon and sample
        // Returns an array of valid abundance values
        let filteredData = [];

        if (sample) {
            // Single sample case
            const sampleData = abundanceDataBySample[sample];
            if (!sampleData || sampleData.length === 0) {
                return filteredData;
            }

            for (const row of sampleData) {
                const matchesTaxon = !taxon || row.taxon === taxon;
                if (matchesTaxon) {
                    const val = typeof row.abundance === 'number' ? row.abundance : parseFloat(row.abundance);
                    if (!isNaN(val) && val > 0) {
                        filteredData.push(val);
                    }
                }
            }
        } else {
            // All samples - aggregate across all loaded samples
            for (const sampleId in abundanceDataBySample) {
                const sampleData = abundanceDataBySample[sampleId];
                for (const row of sampleData) {
                    const matchesTaxon = !taxon || row.taxon === taxon;
                    if (matchesTaxon) {
                        const val = typeof row.abundance === 'number' ? row.abundance : parseFloat(row.abundance);
                        if (!isNaN(val) && val > 0) {
                            filteredData.push(val);
                        }
                    }
                }
            }
        }

        return filteredData;
    }

    function getMeanAbundance(taxon, sample = '') {
        // Use pre-computed mean abundance data from collapsed table
        if (meanAbundancesData && meanAbundancesData[taxon]) {
            if (sample) {
                // Return mean for specific sample
                const value = meanAbundancesData[taxon][sample];
                if (value !== undefined && value !== null) {
                    return value;
                } else {
                    console.error(`Precomputed mean abundance not available for taxon "${taxon}" and sample "${sample}"`);
                    return 0;
                }
            } else {
                // Return overall mean across all samples
                const values = Object.values(meanAbundancesData[taxon]);
                if (values.length > 0) {
                    return d3.mean(values);
                } else {
                    console.error(`Precomputed mean abundance not available for taxon "${taxon}"`);
                    return 0;
                }
            }
        }

        console.error(`Precomputed mean abundance data not available for taxon "${taxon}"`);
        return 0;
    }

    function getContigCount(taxon, sample = '') {
        // Calculate contig count from raw abundance data
        const filteredData = getFilteredAbundanceData(taxon, sample);
        return filteredData.length;
    }

    function hasDataForTaxon(taxon, sampleFilter = '') {
        // Check if there's any data for this taxon (optionally filtered by sample)
        // First try using pre-computed mean abundances (more efficient)
        if (meanAbundancesData && meanAbundancesData[taxon]) {
            if (!sampleFilter) {
                // If no sample filter, any sample with data counts
                return Object.keys(meanAbundancesData[taxon]).length > 0;
            } else {
                // Check if this taxon has data for the specific sample
                return meanAbundancesData[taxon][sampleFilter] !== undefined;
            }
        }

        // Fall back to checking abundanceDataBySample if meanAbundancesData not available
        if (sampleFilter) {
            const sampleData = abundanceDataBySample[sampleFilter];
            if (!sampleData || sampleData.length === 0) {
                return false;
            }
            for (const row of sampleData) {
                if (row.taxon === taxon) {
                    return true;
                }
            }
            return false;
        } else {
            // Check across all samples
            for (const sampleId in abundanceDataBySample) {
                const sampleData = abundanceDataBySample[sampleId];
                for (const row of sampleData) {
                    if (row.taxon === taxon) {
                        return true;
                    }
                }
            }
            return false;
        }
    }

    function calculateStatisticsForSelection(taxon, sample) {
        // Get filtered data for the selected taxon and sample
        const filteredData = getFilteredAbundanceData(taxon, sample);

        if (filteredData.length === 0) {
            return null;
        }

        const count = filteredData.length;
        
        // Use precomputed mean abundance (matches sorting logic)
        let mean;
        if (meanAbundancesData && meanAbundancesData[taxon] && sample) {
            const precomputedMean = meanAbundancesData[taxon][sample];
            if (precomputedMean !== undefined && precomputedMean !== null) {
                mean = precomputedMean;
            } else {
                console.error(`Precomputed mean abundance not available for taxon "${taxon}" and sample "${sample}"`);
                mean = 0;
            }
        } else {
            console.error(`Precomputed mean abundance data not available for taxon "${taxon}"`);
            mean = 0;
        }
        
        const median = d3.median(filteredData);
        const std = d3.deviation(filteredData) || 0;

        return { count, mean, median, std };
    }

    function sortTaxaByValue(taxaToSort, selectedSample, sortBy) {
        // Helper function to filter and sort taxa by a metric (abundance or count)
        // Returns sorted array of taxa
        const taxaWithValues = taxaToSort
            .filter(taxon => hasDataForTaxon(taxon, selectedSample))
            .map(taxon => {
                const value = sortBy === 'count'
                    ? getContigCount(taxon, selectedSample)
                    : getMeanAbundance(taxon, selectedSample);
                return { taxon, value };
            });

        // Sort by selected metric (highest to lowest)
        // Use stable sort: if values are equal, sort by taxon name alphabetically
        taxaWithValues.sort((a, b) => {
            if (b.value !== a.value) {
                return b.value - a.value; // Higher values first
            }
            // If values are equal, sort alphabetically by taxon name
            return a.taxon.localeCompare(b.taxon);
        });
        return taxaWithValues.map(item => item.taxon);
    }

    // ============================================
    // 5. UI UPDATE FUNCTIONS
    // ============================================

    function updateSelectionStats() {
        const selectedTaxonLabel = taxonDropdown.value;
        const selectedSample = sampleDropdown.value || '';

        // If no taxon is selected, show default message
        if (!selectedTaxonLabel) {
            selectionStats.innerHTML = `
                <small class="text-muted">Select taxon and sample to view statistics</small>
            `;
            return;
        }

        const selectedTaxon = taxaMap[selectedTaxonLabel] || selectedTaxonLabel;
        const stats = calculateStatisticsForSelection(selectedTaxon, selectedSample);

        if (!stats) {
            selectionStats.innerHTML = `
                <small class="text-muted">No data available for this selection</small>
            `;
            return;
        }

        selectionStats.innerHTML = `
            <div class="stats-content">
                <span class="stat-item">Contig count: <span class="badge bg-success">${stats.count.toLocaleString()}</span></span>
                <span class="stat-item">Mean abundance: <span class="badge bg-success">${stats.mean.toFixed(1)}</span></span>
                <span class="stat-item">Median abundance: <span class="badge bg-success">${stats.median.toFixed(1)}</span></span>
                <span class="stat-item">Abundance std dev: <span class="badge bg-success">${stats.std.toFixed(1)}</span></span>
            </div>
        `;
    }

    function updateTaxonDropdown(selectedSample) {
        // Filter and populate taxon dropdown based on selected sample
        if (!meanAbundancesData || !taxaList || taxaList.length === 0) return;

        // Get currently selected taxon to preserve selection if possible
        const currentSelection = taxonDropdown.value;

        // Clear dropdown and add "None" option
        taxonDropdown.innerHTML = '<option value="">None</option>';

        // Filter taxa that have data for the selected sample
        const availableTaxonLabels = [];
        taxaMap = {}; // Reset map

        taxaList.forEach(taxon => {
            // Check if taxon has data for selected sample
            if (hasDataForTaxon(taxon, selectedSample)) {
                const taxonLabel = getTaxonDropdownLabel(taxon);
                taxaMap[taxonLabel] = taxon;
                availableTaxonLabels.push(taxonLabel);
            }
        });

        // Sort taxon labels alphabetically for dropdown
        availableTaxonLabels.sort();

        // Populate dropdown
        availableTaxonLabels.forEach(taxonLabel => {
            const opt = document.createElement('option');
            opt.value = taxonLabel;
            opt.textContent = taxonLabel;
            taxonDropdown.appendChild(opt);
        });

        // Try to restore previous selection if it's still available
        if (currentSelection && availableTaxonLabels.includes(currentSelection)) {
            taxonDropdown.value = currentSelection;
        } else {
            taxonDropdown.value = '';
        }
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

    function resortTaxaByMeanAbundance(selectedSample) {
        // Resort taxa by mean abundance or contig count based on user selection
        if (!meanAbundancesData || taxaList.length === 0) return;

        const sortBy = document.querySelector('input[name="sortBy"]:checked')?.value || 'abundance';

        // Sort taxa using shared helper
        const newTaxaList = sortTaxaByValue(taxaList, selectedSample, sortBy);

        // Reorder DOM elements to match new sort order
        const histogramItems = Array.from(document.querySelectorAll('.histogram-item'));
        const itemsByTaxon = {};
        histogramItems.forEach(item => {
            itemsByTaxon[item.dataset.taxon] = item;
        });

        // Create a fragment and append items in new order
        // appendChild automatically moves elements if they're already in the DOM
        const fragment = document.createDocumentFragment();
        newTaxaList.forEach(taxon => {
            if (itemsByTaxon[taxon]) {
                fragment.appendChild(itemsByTaxon[taxon]);
            }
        });

        // Clear grid and append fragment (this preserves plot containers)
        histogramGrid.innerHTML = '';
        histogramGrid.appendChild(fragment);

        taxaList = newTaxaList;
    }

    // ============================================
    // 6. PLOT CREATION & UPDATE FUNCTIONS
    // ============================================

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

        // Create a spinner for loading indication
        const spinnerDiv = document.createElement('div');
        spinnerDiv.className = 'histogram-spinner';
        spinnerDiv.id = `spinner-${taxon.replace(/[^a-zA-Z0-9]/g, '_')}`;
        spinnerDiv.innerHTML = '<div class="rect1"></div><div class="rect2"></div><div class="rect3"></div><div class="rect4"></div><div class="rect5"></div>';
        containerDiv.appendChild(spinnerDiv);

        // Create a div for the Vega visualization or "no data" message
        // Don't create the plot yet - it will be created lazily in updateVisualization()
        const plotDiv = document.createElement('div');
        plotDiv.id = `histogram-${taxon.replace(/[^a-zA-Z0-9]/g, '_')}`;
        plotDiv.className = 'histogram-plot';
        containerDiv.appendChild(plotDiv);

        histogramGrid.appendChild(containerDiv);

        // Don't create plots here - let updateVisualization() handle it lazily
        // This ensures we only create plots for visible items
    }

    function createVegaPlot(taxon, plotDiv, spinner, selectedSample) {
        // Helper function to create a Vega plot for a given taxon
        plotDiv.className = 'histogram-plot';
        plotDiv.innerHTML = '';

        // Show spinner
        if (spinner) {
            spinner.style.display = 'block';
        }

        // Load JSON data for the selected sample
        return loadSampleData(selectedSample).then(() => {
            const sampleData = abundanceDataBySample[selectedSample] || [];

            const spec = JSON.parse(JSON.stringify(vegaAbundanceHistogramSpec));
            spec.params[0].value = taxon;

            return vegaEmbed(plotDiv, spec, {
                actions: false,
                renderer: 'canvas'
            }).then(res => {
                vegaViews[taxon] = res.view;
                // Set the data and run
                return res.view.data('source', sampleData).runAsync().then(() => {
                    // Hide spinner when plot is loaded
                    if (spinner) {
                        spinner.style.display = 'none';
                    }
                    // Force a resize after a short delay to ensure grid layout is complete
                    setTimeout(() => {
                        try {
                            res.view.resize();
                            res.view.runAsync();
                        } catch (e) {
                            console.warn(`Could not resize view for ${taxon}:`, e);
                        }
                    }, 100);
                });
            }).catch(error => {
                console.error(`Error embedding histogram for taxon ${taxon}:`, error);
                plotDiv.innerHTML = '<p class="text-danger" style="font-size: 10px;">Error loading plot</p>';
                // Hide spinner on error
                if (spinner) {
                    spinner.style.display = 'none';
                }
                // Clean up view reference on error
                delete vegaViews[taxon];
            });
        }).catch(error => {
            console.error(`Error loading sample data for ${selectedSample}:`, error);
            plotDiv.innerHTML = '<p class="text-danger" style="font-size: 10px;">Error loading data</p>';
            if (spinner) {
                spinner.style.display = 'none';
            }
        });
    }

    function updateVisualization() {
        const selectedSample = sampleDropdown.value || '';
        const taxaLimit = taxaLimitSlider ? parseInt(taxaLimitSlider.value) : taxaList.length;

        // Create a map of taxon to DOM element for efficient lookup
        const itemsByTaxon = {};
        const histogramItems = document.querySelectorAll('.histogram-item');
        histogramItems.forEach(item => {
            const taxon = item.dataset.taxon;
            if (taxon) {
                itemsByTaxon[taxon] = item;
            }
        });

        // Iterate through taxaList (the sorted order) instead of DOM order
        taxaList.forEach((taxon, index) => {
            const item = itemsByTaxon[taxon];
            if (!item) return;

            const plotDiv = item.querySelector('.histogram-plot, .histogram-no-data');
            if (!plotDiv) return;

            // Hide items beyond the limit
            if (index >= taxaLimit) {
                item.style.display = 'none';
                // Destroy view for hidden plots to free up memory
                if (vegaViews[taxon]) {
                    try {
                        vegaViews[taxon].finalize();
                    } catch (e) {
                        // Ignore errors during finalization
                    }
                    delete vegaViews[taxon];
                }
                return;
            }

            // Check if there's data for this taxon with the current sample filter
            if (!hasDataForTaxon(taxon, selectedSample)) {
                // Hide taxa with no data (no checkbox option)
                item.style.display = 'none';
                // Destroy view if it exists
                if (vegaViews[taxon]) {
                    try {
                        vegaViews[taxon].finalize();
                    } catch (e) {
                        // Ignore errors during finalization
                    }
                    delete vegaViews[taxon];
                }
                return;
            } else {
                // There is data - show plot
                item.style.display = '';

                // Get spinner element for this taxon
                const spinnerId = `spinner-${taxon.replace(/[^a-zA-Z0-9]/g, '_')}`;
                const spinner = item.querySelector(`#${spinnerId}`);

                // Check if plot already exists and is still valid
                const plotExists = plotDiv.querySelector('canvas, svg');

                if (!vegaViews[taxon] || !plotExists || plotDiv.className === 'histogram-no-data') {
                    // Plot doesn't exist yet, was destroyed, or needs recreation - create it
                    createVegaPlot(taxon, plotDiv, spinner, selectedSample);
                } else {
                    // Update existing plot with new sample data (only if visible)
                    // Check if view is still valid before updating
                    if (vegaViews[taxon] && plotDiv.querySelector('canvas, svg')) {
                        loadSampleData(selectedSample).then(() => {
                            const sampleData = abundanceDataBySample[selectedSample] || [];
                            try {
                                vegaViews[taxon].data('source', sampleData).runAsync();
                            } catch (e) {
                                // View might have been destroyed, recreate it
                                console.warn(`View for ${taxon} was invalid, recreating...`);
                                delete vegaViews[taxon];
                                // Trigger recreation by clearing plotDiv
                                plotDiv.innerHTML = '';
                                // This will be handled in the next updateVisualization call
                            }
                        }).catch(error => {
                            console.error(`Error loading sample data for ${selectedSample}:`, error);
                        });
                    }
                }
            }
        });

        // Hide any items in DOM that are not in taxaList (shouldn't happen, but safety check)
        histogramItems.forEach(item => {
            const taxon = item.dataset.taxon;
            if (taxon && !taxaList.includes(taxon)) {
                item.style.display = 'none';
            }
        });
    }

    // ============================================
    // 7. EVENT LISTENERS
    // ============================================

    // Event listener for sample dropdown
    sampleDropdown.addEventListener('change', () => {
        const selectedSample = sampleDropdown.value || '';
        // Load the new sample data first, then update everything
        loadSampleData(selectedSample).then(() => {
            updateTaxonDropdown(selectedSample);
            resortTaxaByMeanAbundance(selectedSample);
            debouncedUpdateVisualization();
            updateSelectionStats();
        }).catch(error => {
            console.error(`Error loading sample data for ${selectedSample}:`, error);
        });
    });

    // Event listener for taxon dropdown
    taxonDropdown.addEventListener('change', updateTaxonHighlight);

    // Event listener for taxa limit slider
    if (taxaLimitSlider) {
        taxaLimitSlider.addEventListener('input', () => {
            const value = parseInt(taxaLimitSlider.value);
            if (taxaLimitValue) {
                taxaLimitValue.textContent = value;
            }
            debouncedUpdateVisualization();
        });
    }

    // Event listener for sort by radio buttons
    const sortByRadios = document.querySelectorAll('input[name="sortBy"]');
    sortByRadios.forEach(radio => {
        radio.addEventListener('change', () => {
            const selectedSample = sampleDropdown.value || '';
            resortTaxaByMeanAbundance(selectedSample);
            debouncedUpdateVisualization();
        });
    });

    // ============================================
    // 8. INITIALIZATION
    // ============================================

    // Populate sample dropdown - default to first sample (no "All" option)
    if (samples && Array.isArray(samples) && samples.length > 0) {
        // Sort samples alphabetically
        const sortedSamples = [...samples].sort();
        sortedSamples.forEach((sample, index) => {
            const opt = document.createElement('option');
            opt.value = sample;
            opt.textContent = sample;
            sampleDropdown.appendChild(opt);
            // Set first sample as default
            if (index === 0) {
                sampleDropdown.value = sample;
            }
        });
    } else {
        console.warn('No samples data available for dropdown');
    }


    // Load initial sample data from JSON file
    const initialSample = sampleDropdown.value || (samples && samples.length > 0 ? samples[0] : '');
    if (!initialSample) {
        console.error('No samples available');
        document.getElementById('spinner-histogram')?.remove();
        histogramGrid.innerHTML = '<p class="text-danger">Error: No samples available</p>';
        return;
    }

    loadSampleData(initialSample)
        .then(() => {
            // Verify data is loaded
            if (!abundanceDataBySample[initialSample] || abundanceDataBySample[initialSample].length === 0) {
                console.error('No data loaded for initial sample:', initialSample);
                document.getElementById('spinner-histogram')?.remove();
                histogramGrid.innerHTML = '<p class="text-danger">Error: No data loaded for sample</p>';
                return;
            }

            // Extract unique taxa from the data
            const taxaSet = new Set();
            for (const sampleId in abundanceDataBySample) {
                const sampleData = abundanceDataBySample[sampleId];
                for (const row of sampleData) {
                    if (row.taxon) {
                        taxaSet.add(row.taxon);
                    }
                }
            }

            // Get selected sample for sorting
            const selectedSample = sampleDropdown.value || '';
            const sortBy = document.querySelector('input[name="sortBy"]:checked')?.value || 'abundance';

            // Sort taxa using shared helper
            taxaList = sortTaxaByValue(Array.from(taxaSet), selectedSample, sortBy);

            // Set slider max value if it exists
            if (taxaLimitSlider) {
                taxaLimitSlider.max = taxaList.length;
                taxaLimitSlider.value = Math.min(20, taxaList.length);
                if (taxaLimitValue) {
                    taxaLimitValue.textContent = taxaLimitSlider.value;
                }
            }

            // Populate taxon dropdown based on selected sample
            // (will be filtered to only show taxa with data for selected sample)
            updateTaxonDropdown(selectedSample);

            // Remove spinner
            document.getElementById('spinner-histogram')?.remove();

            // Create a histogram for each taxon (only create DOM elements, plots will be created lazily)
            if (Object.keys(vegaAbundanceHistogramSpec).length > 0 && taxaList.length > 0) {
                taxaList.forEach(taxon => {
                    const taxonShort = getTaxonShort(taxon);
                    createHistogramForTaxon(taxon, taxonShort);
                });

                // Wait for DOM to finish layout before creating plots
                // This ensures containers have correct dimensions
                setTimeout(() => {
                    updateVisualization();
                    updateSelectionStats();
                }, 100);
            } else {
                histogramGrid.innerHTML = '<p class="text-danger">Error: No data or invalid visualization specification.</p>';
            }
        })
        .catch(error => {
            console.error('Error loading initial sample data:', error);
            document.getElementById('spinner-histogram')?.remove();
            histogramGrid.innerHTML = '<p class="text-danger">Error loading data: ' + (error.message || error) + '</p>';
        });
});
