// Plotly configuration
const plotConfig = {
    responsive: true,
    displayModeBar: true,
    displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    toImageButtonOptions: {
        format: 'png',
        filename: 'fep_analysis_plot',
        height: 600,
        width: 1000,
        scale: 2
    }
};

function buildAxisStyle(title) {
    return {
        title: {
            text: `<b>${title}</b>`,
            font: { size: 16, color: '#222' }
        },
        tickfont: { size: 14, color: '#333' },
        gridcolor: '#e0e0e0',
        zeroline: false
    };
}

function buildBaseLayout(title, xTitle, yTitle, extra = {}) {
    return Object.assign({
        title: { text: `<b>${title}</b>`, font: { size: 16, color: '#222' } },
        xaxis: buildAxisStyle(xTitle),
        yaxis: buildAxisStyle(yTitle),
        plot_bgcolor: 'rgba(255,255,255,0.8)',
        paper_bgcolor: 'rgba(0,0,0,0)',
        margin: { l: 78, r: 30, t: 52, b: 62 },
        font: { family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif', size: 14, color: '#333' },
        autosize: true
    }, extra);
}

// Initialize on DOM ready
document.addEventListener('DOMContentLoaded', function() {
    // Animate cards on scroll
    const observerOptions = {
        threshold: 0.1,
        rootMargin: '0px 0px -50px 0px'
    };

    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.style.animationPlayState = 'running';
            }
        });
    }, observerOptions);

    document.querySelectorAll('.card').forEach(card => {
        observer.observe(card);
    });

    // Add click-to-copy for values
    document.querySelectorAll('.summary-item').forEach(item => {
        item.addEventListener('click', function() {
            const value = this.querySelector('.value');
            if (value) {
                const text = value.textContent.trim();
                navigator.clipboard.writeText(text).then(() => {
                    // Show feedback
                    const feedback = document.createElement('div');
                    feedback.textContent = '✓ Copied!';
                    feedback.style.cssText = 'position:fixed;top:20px;right:20px;background:#4CAF50;color:white;padding:10px 20px;border-radius:5px;animation:fadeInOut 2s forwards;';
                    document.body.appendChild(feedback);
                    setTimeout(() => feedback.remove(), 2000);
                });
            }
        });
    });
});

// Create Plotly plot
function createPlot(divId, data, layout) {
    Plotly.newPlot(divId, data, layout, plotConfig);
}

// Tab switching
function showBackend(backendId) {
    // Hide all tabs
    document.querySelectorAll('.tab-content').forEach(el => el.style.display = 'none');
    document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));

    // Show selected tab
    const targetTab = document.getElementById('backend-' + backendId);
    if (targetTab) {
        targetTab.style.display = 'block';
    }

    // Update active button
    const activeBtn = document.querySelector(`[onclick="showBackend('${backendId}')"]`);
    if (activeBtn) {
        activeBtn.classList.add('active');
    }
}

// Bound/Unbound toggle
function updateLambdaPlots(backendId) {
    const showBound = document.getElementById('show-bound-' + backendId)?.checked ?? true;
    const showUnbound = document.getElementById('show-unbound-' + backendId)?.checked ?? true;

    // Update Plotly traces visibility
    const plotDiv = document.getElementById('dg-lambda-plot-' + backendId);
    if (plotDiv && plotDiv.data) {
        plotDiv.data.forEach((trace, index) => {
            const isBound = trace.name.includes('Bound');
            const isVisible = (isBound && showBound) || (!isBound && showUnbound);
            Plotly.restyle(plotDiv, {visible: isVisible}, [index]);
        });
    }
}

// Generate comparison plot
function generateComparisonPlot(containerId, data) {
    const trace = {
        x: data.names,
        y: data.values,
        error_y: {
            type: 'data',
            array: data.errors,
            visible: true,
            color: 'rgba(100, 100, 100, 0.5)',
            thickness: 1.5,
            width: 10
        },
        type: 'bar',
        marker: {
            color: ['rgb(76, 175, 80)', 'rgb(33, 150, 243)', 'rgb(255, 152, 0)']
        },
        text: data.values.map(v => v.toFixed(2)),
        textposition: 'auto',
        textfont: { size: 12, color: '#333' }
    };

    const layout = buildBaseLayout('Binding Free Energy (ΔG) Comparison', 'Backend / Estimator', 'ΔG (kcal/mol)', {
        margin: { l: 78, r: 30, t: 52, b: 100 },
        showlegend: false
    });
    layout.xaxis.tickangle = -45;

    Plotly.newPlot(containerId, [trace], layout, plotConfig);
}

// Generate overlap matrix heatmap
function generateOverlapMatrix(containerId, matrix, lambdaValues) {
    const data = [{
        z: matrix,
        x: lambdaValues.map(v => v.toFixed(2)),
        y: lambdaValues.map(v => v.toFixed(2)),
        type: 'heatmap',
        colorscale: 'RdYlGn',
        reversescale: false,
        colorbar: {
            title: {
                text: 'Overlap',
                font: { size: 14 }
            }
        },
        hovertemplate: '<b>λ %{x}</b> → <b>λ %{y}</b><br>Overlap: %{z:.3f}<extra></extra>'
    }];

    const layout = buildBaseLayout('Lambda State Overlap Matrix', 'Lambda λ', 'Lambda λ', {
        margin: { l: 78, r: 60, t: 52, b: 62 },
        height: 420,
        hovermode: 'closest'
    });

    Plotly.newPlot(containerId, data, layout, plotConfig);
}

// Estimator tab switching for multi-estimator reports
function showEstimatorTab(estimatorId) {
    'use strict';
    // Hide all tab contents
    const tabContents = document.querySelectorAll('.tab-content');
    tabContents.forEach(el => {
        el.style.display = 'none';
    });

    // Remove active class from all buttons
    const tabBtns = document.querySelectorAll('.tab-btn');
    tabBtns.forEach(btn => {
        btn.classList.remove('active');
    });

    // Show selected tab content
    const targetTab = document.getElementById('estimator-' + estimatorId);
    if (targetTab) {
        targetTab.style.display = 'block';
    }

    // Add active class to clicked button
    const activeBtn = document.querySelector(`[onclick="showEstimatorTab('${estimatorId}')"]`);
    if (activeBtn) {
        activeBtn.classList.add('active');
    }

    // Refresh Plotly charts in visible tab (if any)
    // This ensures charts resize properly when tab becomes visible
    if (targetTab && window.Plotly) {
        window.setTimeout(() => {
            const plotDivs = targetTab.querySelectorAll('.js-plotly-plot, [id*="plot"], [id*="time-conv"], [id*="bootstrap-hist"], [id*="overlap-matrix"]');
            plotDivs.forEach(div => {
                try {
                    Plotly.Plots.resize(div);
                    Plotly.relayout(div, { autosize: true });
                } catch (err) {
                    // Ignore non-plot containers.
                }
            });
        }, 80);
    }
}
