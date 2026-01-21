// Beam selector event listener
const nbeamsSelector = document.querySelector("#beamSelector");
nbeamsSelector.addEventListener("change", function () {
    const inputs = document.getElementsByClassName("type1");
    const nbeams = parseInt(nbeamsSelector.value);
    for (let i = 0; i < inputs.length; i++) {
        inputs[i].style.display = i < nbeams ? 'grid' : 'none';
    }
    mainFunction();
});

// Point loads selector event listener
const ploadsSelector = document.getElementById("ploadSelector");
ploadsSelector.addEventListener("change", function () {
    const inputs = document.getElementsByClassName("type2");
    const npoints = parseInt(ploadsSelector.value);
    for (let i = 0; i < inputs.length; i++) {
        inputs[i].style.display = i < npoints ? 'grid' : 'none';
    }
    mainFunction();
});

// Attach change listeners to all inputs
$(".lengths, .loads, .ploads, .ploadsx, .custom-select").change(mainFunction);

// Run main function on startup
mainFunction();

function mainFunction() {
    // Get visible counts from selectors
    const nbeams = parseInt(nbeamsSelector.value);
    const npoints = parseInt(ploadsSelector.value);

    // Collect beam lengths and loads using IDs
    const L = [];
    const loads = [];
    for (let i = 1; i <= nbeams; i++) {
        L.push(parseFloat(document.getElementById('L' + i).value));
        loads.push(parseFloat(document.getElementById('q' + i).value) / 100);
    }

    // Calculate node positions
    const nodesX = [0];
    for (let i = 0; i < L.length; i++) {
        nodesX[i + 1] = nodesX[i] + L[i];
    }

    // Update max value for point load positions
    const totalLength = nodesX[nodesX.length - 1];
    for (let i = 1; i <= 3; i++) {
        document.getElementById('x' + i).max = totalLength;
    }

    // Collect point loads using IDs
    const pointLoads = { x: [], y: [] };
    for (let i = 1; i <= npoints; i++) {
        pointLoads.y.push(parseFloat(document.getElementById('P' + i).value) / 100);
        pointLoads.x.push(parseFloat(document.getElementById('x' + i).value));
    }

    // Get boundary conditions
    const BC = new Array(nodesX.length).fill(1);
    BC[0] = parseInt(document.getElementById('LeftSupportCondition').value);
    BC[BC.length - 1] = parseInt(document.getElementById('RightSupportCondition').value);

    // Calculate segment lengths and total length
    const segmentLengths = [];
    for (let j = 0; j < nodesX.length - 1; j++) {
        segmentLengths[j] = nodesX[j + 1] - nodesX[j];
    }
    const totLength = segmentLengths.reduce((curr, prev) => curr + prev, 0);

    // Determine constrained DOFs and zero reactions
    const u_zero = [];
    const R_zero = [];

    for (let i = 0; i < BC.length; i++) {
        if (BC[i] === 0) {
            // Free: both reactions are zero
            R_zero.push(i * 2);
            R_zero.push(i * 2 + 1);
        } else if (BC[i] === 1) {
            // Pinned: vertical displacement is zero, moment reaction is zero
            u_zero.push(i * 2);
            R_zero.push(i * 2 + 1);
        } else if (BC[i] === 2) {
            // Rigid: both displacements are zero
            u_zero.push(i * 2);
            u_zero.push(i * 2 + 1);
        }
    }

    // Draw geometry first
    geoDraw(nodesX, loads, BC, pointLoads, totLength);

    // Build element stiffness matrices
    const KeAll = [];
    for (let i = 0; i < nodesX.length - 1; i++) {
        const Li = segmentLengths[i];
        const Li3 = Li * Li * Li;
        const Li2 = Li * Li;

        KeAll[i] = [
            [12 / Li3, 6 * Li / Li3, -12 / Li3, 6 * Li / Li3],
            [6 * Li / Li3, 4 * Li2 / Li3, -6 * Li / Li3, 2 * Li2 / Li3],
            [-12 / Li3, -6 * Li / Li3, 12 / Li3, -6 * Li / Li3],
            [6 * Li / Li3, 2 * Li2 / Li3, -6 * Li / Li3, 4 * Li2 / Li3]
        ];
    }

    // Assemble global stiffness matrix
    function compile(KeAll) {
        const sz = 2 + KeAll.length * 2;
        const Ks = Array.from({ length: sz }, () => new Array(sz).fill(0));

        for (let i = 0; i < KeAll.length; i++) {
            const idx = i * 2;
            for (let j = 0; j < 4; j++) {
                for (let k = 0; k < 4; k++) {
                    Ks[idx + j][idx + k] += KeAll[i][j][k];
                }
            }
        }
        return Ks;
    }

    let Ks = compile(KeAll);
    const sz = Ks.length;

    // Round values for numerical stability
    for (let i = 0; i < sz; i++) {
        Ks[i] = Ks[i].map(x => Math.round(x * 10000000) / 10000000);
    }

    // Create symbolic variables
    const u = [];
    const R = [];
    for (let i = 0; i < sz; i++) {
        u[i] = 'x' + (i + 1);
        R[i] = 'y' + (i + 1);
    }

    // Calculate equivalent nodal forces from distributed loads
    const F_n = [];
    for (let i = 0; i < nodesX.length; i++) {
        if (i === 0) {
            F_n[2 * i] = -loads[i] * segmentLengths[0] / 2;
            F_n[2 * i + 1] = -loads[i] * segmentLengths[0] ** 2 / 12;
        } else if (i === nodesX.length - 1) {
            F_n[i * 2] = -loads[i - 1] * segmentLengths[i - 1] / 2;
            F_n[i * 2 + 1] = loads[i - 1] * segmentLengths[i - 1] ** 2 / 12;
        } else {
            F_n[i * 2] = -loads[i - 1] * segmentLengths[i - 1] / 2 - loads[i] * segmentLengths[i] / 2;
            F_n[i * 2 + 1] = loads[i - 1] * segmentLengths[i - 1] ** 2 / 12 - loads[i] * segmentLengths[i] ** 2 / 12;
        }
    }

    // Add point load contributions to nodal forces
    for (let j = 0; j < pointLoads.x.length; j++) {
        const id_end = nodesX.map((x, i) => x >= pointLoads.x[j] ? i : null).filter(x => x !== null)[0];
        const id_start = id_end - 1;

        const x_1 = pointLoads.x[j] - nodesX[id_start];
        const x_2 = nodesX[id_end] - nodesX[id_start] - x_1;
        const Li = x_1 + x_2;
        const Li2 = Li * Li;
        const Li3 = Li * Li * Li;

        F_n[2 * id_start] -= pointLoads.y[j] * (3 * x_1 + x_2) * x_2 * x_2 / Li3;
        F_n[2 * id_start + 1] -= pointLoads.y[j] * x_2 * x_2 * x_1 / Li2;
        F_n[2 * id_end] -= pointLoads.y[j] * (3 * x_2 + x_1) * x_1 * x_1 / Li3;
        F_n[2 * id_end + 1] += pointLoads.y[j] * x_1 * x_1 * x_2 / Li2;
    }

    // Build system of equations
    const eq = [];
    for (let j = 0; j < sz; j++) {
        let eqStr = R[j] + '+' + F_n[j] + '=';
        for (let i = 0; i < sz; i++) {
            eqStr += '+' + Ks[j][i] + '*' + u[i];
        }
        eq.push(eqStr);
    }

    // Add boundary conditions
    for (let i = 0; i < u_zero.length; i++) {
        eq.push(u[u_zero[i]] + '=0');
    }
    for (let i = 0; i < R_zero.length; i++) {
        eq.push(R[R_zero[i]] + '=0');
    }

    // Solve system using Nerdamer
    const sol2 = nerdamer.solveEquations(eq, [...u, ...R]);
    console.log(sol2.toString());

    // Helper function for linspace
    function linspace(min, max, n) {
        const nel = n - 1;
        const step = (max - min) / nel;
        return Array.from({ length: n }, (_, i) => min + i * step);
    }

    // Draw geometry
    function geoDraw(nodesX, loads, BC, pointLoads, totLength) {
        // Scale factor for supports - keeps them visually constant size
        const s = totLength * 0.02; // Support size scales with total length
        const margin = totLength * 0.02; // Margin for axis range

        const path = [];
        for (let j = 0; j < nodesX.length; j++) {
            const x = nodesX[j];
            if (BC[j] === 1) {
                // Pinned (triangle)
                path[j] = 'M' + x + ',0L' + (x + s) + ',-' + s + 'L' + (x - s) + ',-' + s + ' Z';
            } else if (BC[j] === 2) {
                // Rigid (rectangle)
                if (j + 1 < nodesX.length) {
                    path[j] = 'M' + x + ',' + s + 'L' + x + ',-' + s + 'L' + (x - s) + ',-' + s + 'L' + (x - s) + ',' + s + ' Z';
                } else {
                    path[j] = 'M' + x + ',' + s + 'L' + x + ',-' + s + 'L' + (x + s) + ',-' + s + 'L' + (x + s) + ',' + s + ' Z';
                }
            } else if (BC[j] === 0) {
                // Free
                path[j] = 'M' + x + ',0L0,0 Z';
            }
        }

        const shapes = [];
        for (let j = 0; j < nodesX.length; j++) {
            shapes[j] = {
                type: 'path',
                path: path[j],
                fillcolor: 'rgba(148, 163, 184, 0.5)',
                line: { color: 'rgb(71, 85, 105)', width: 1 }
            };
        }

        // UDL load rectangles - scale height with total length
        const loadShapes = [];
        for (let i = 0; i < nodesX.length - 1; i++) {
            loadShapes[i] = {
                type: 'rect',
                x0: nodesX[i],
                x1: nodesX[i + 1],
                y0: s,
                y1: loads[i] * totLength * 0.1 + s,
                fillcolor: 'rgba(34, 197, 94, 0.4)',
                line: { width: 0 }
            };
        }

        // Point load annotations
        const ploadAnnotations = [];
        for (let i = 0; i < pointLoads.x.length; i++) {
            ploadAnnotations[i] = {
                x: pointLoads.x[i],
                y: 0.01,
                xref: 'x',
                yref: 'y',
                text: 'P<sub>' + (i + 1) + '</sub>',
                font: { size: 20, color: '#dc2626' },
                arrowcolor: '#dc2626',
                showarrow: true,
                arrowhead: 9,
                ax: 0,
                ay: -50 - (pointLoads.y[i] * 100)
            };
        }

        const allShapes = [...shapes, ...loadShapes];

        const layout = {
            annotations: ploadAnnotations,
            shapes: allShapes,
            xaxis: {
                scaleanchor: "y",
                range: [-margin, totLength + margin],
                fixedrange: true,
                tickvals: nodesX,
                zeroline: false,
                gridcolor: 'rgba(148, 163, 184, 0.2)'
            },
            yaxis: {
                visible: true,
                showgrid: false,
                fixedrange: true
            },
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            margin: { l: 0, b: 0, t: 10, r: 0 },
            height: 280
        };

        const trace1 = {
            x: nodesX,
            y: new Array(nodesX.length).fill(0),
            line: { color: '#1e293b', width: 4 },
            marker: { size: 8, color: '#1e293b' },
            hoverinfo: 'x'
        };

        const myAx = document.getElementById("myFig");
        Plotly.newPlot(myAx, [trace1], layout, { responsive: true, displayModeBar: false });
    }

    // Calculate shear force, moment, and deflection diagrams
    const x_tot = [];
    const y_tot = [];
    const M_tot = [];
    const v_tot = [];
    const v_right = [sol2[1][1]];

    // Calculate cumulative reaction forces
    const Ri = [sol2[sz][1]];
    for (let i = 1; i < sz / 2; i++) {
        Ri[i] = Ri[i - 1] + sol2[sz + i * 2][1];
    }

    // Process each beam segment
    const Vi = [];
    for (let i = 0; i < nodesX.length - 1; i++) {
        const x_start = nodesX[i];
        const x_end = nodesX[i + 1];
        const n_el = 10000;
        const x = linspace(x_start, x_end, n_el);
        const dx = (x_end - x_start) / n_el;

        // Find point loads to the left of current segment
        const pID_left = pointLoads.x.map((px, idx) => px < x_start ? idx : null).filter(idx => idx !== null);

        let pSumLeft = 0;
        if (pID_left.length > 0) {
            pSumLeft = pID_left.map(idx => pointLoads.y[idx]).reduce((curr, prev) => curr + prev);
        }

        // Sum of distributed loads to the left
        let qSumLeft = 0;
        for (let k = 0; k < i; k++) {
            qSumLeft += loads[k] * (nodesX[k + 1] - nodesX[k]);
        }

        Vi[i] = Ri[i] - qSumLeft - pSumLeft;

        // Initialize shear force array
        let y = new Array(x.length).fill(Vi[i]);

        // Find point loads within current segment
        const pID_seg = pointLoads.x.map((px, idx) => (px >= x_start && px < x_end) ? idx : null).filter(idx => idx !== null);

        // Apply point load discontinuities
        for (let m = 0; m < pID_seg.length; m++) {
            y = x.map((xi, k) => xi >= pointLoads.x[pID_seg[m]] ? y[k] - pointLoads.y[pID_seg[m]] : y[k]);
        }

        // Apply distributed load
        y = y.map((yi, k) => yi - loads[i] * dx * k);

        // Integrate shear to get moment
        const M = new Array(x.length);
        if (i === 0) {
            M[0] = -sol2[sz + 1][1]; // moment at left reaction
        } else {
            M[0] = M_tot[M_tot.length - 1]; // moment from left side of support
        }
        for (let j = 1; j < x.length; j++) {
            M[j] = M[j - 1] + dx * y[j];
        }

        // Double integrate moment to get deflection
        const a = new Array(x.length);
        const v = new Array(x.length);
        a[0] = sol2[i * 2 + 1][1];
        v[0] = sol2[i * 2][1];
        for (let j = 1; j < x.length; j++) {
            a[j] = a[j - 1] + dx * M[j];
            v[j] = v[j - 1] + dx * a[j];
        }

        x_tot.push(...x);
        y_tot.push(...y);
        M_tot.push(...M);
        v_tot.push(...v);
        v_right.push(v[v.length - 1]);
    }

    // Plot results
    const ax2 = document.getElementById('myPlot2');

    // Calculate scale factors
    const k_M = Math.max(Math.max(...M_tot), Math.abs(Math.min(...M_tot))) || 1;
    const k_V = Math.max(Math.max(...y_tot), Math.abs(Math.min(...y_tot))) || 1;
    const k_v = Math.max(Math.max(...v_tot), Math.abs(Math.min(...v_tot))) || 1;

    // Check for numerical errors
    const error = v_right.map((v, i) => Math.abs((v - sol2[i * 2][1])) / k_v).some(x => x > 0.05);

    const trace_V = {
        x: x_tot,
        y: y_tot.map(x => x / k_V),
        mode: 'lines',
        fill: 'tozeroy',
        fillcolor: 'rgba(59, 130, 246, 0.2)',
        line: { color: '#3b82f6', width: 2 },
        name: 'Shear force, scale: ' + k_V.toExponential(2)
    };

    const trace_M = {
        x: x_tot,
        y: M_tot.map(x => x / k_M),
        mode: 'lines',
        fill: 'tozeroy',
        fillcolor: 'rgba(16, 185, 129, 0.2)',
        line: { color: '#10b981', width: 2 },
        name: 'Moment, scale: ' + k_M.toExponential(2)
    };

    const trace_v = {
        x: x_tot,
        y: v_tot.map(x => x / k_v),
        mode: 'lines',
        line: { color: '#f59e0b', width: 2 },
        name: 'Deflection, scale: ' + k_v.toExponential(2)
    };

    const annotations = [{
        xref: 'paper',
        yref: 'paper',
        x: 0.5,
        xanchor: 'center',
        y: 0.9,
        yanchor: 'bottom',
        text: 'ups... large error detected',
        showarrow: false,
        font: { color: '#dc2626', size: 24 },
        opacity: error ? 1 : 0
    }];

    const layout2 = {
        xaxis: {
            range: [-0.1, totLength + 0.1],
            title: "",
            fixedrange: true,
            tickvals: nodesX,
            zeroline: false,
            gridcolor: 'rgba(148, 163, 184, 0.2)'
        },
        yaxis: {
            title: "c(x)",
            fixedrange: true,
            showgrid: false,
            zeroline: true,
            zerolinecolor: 'rgba(148, 163, 184, 0.5)'
        },
        height: 220,
        showlegend: true,
        legend: {
            orientation: "h",
            y: -0.2,
            font: { family: 'system-ui, sans-serif', size: 11, color: '#64748b' },
            bgcolor: 'rgba(255,255,255,0.9)',
            bordercolor: '#e2e8f0',
            borderwidth: 1
        },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        margin: { l: 0, b: 60, t: 10, r: 0 },
        annotations: annotations,
        autosize: true
    };

    const config = { responsive: true, displayModeBar: false };
    const data2 = [trace_M, trace_V, trace_v];
    Plotly.newPlot(ax2, data2, layout2, config);

    // Update stiffness matrix display
    const p = document.getElementById('modal-body-text');
    let matrixString = '$$EI\\begin{bmatrix}';

    for (let i = 0; i < Ks.length; i++) {
        const rowString = Ks[i].map(x => Math.round(x)).toString().replaceAll(',', '&');
        if (i < Ks.length) {
            matrixString = matrixString.concat(rowString, '\\\\');
        }
    }

    matrixString = matrixString.concat('\\end{bmatrix}', '\\begin{bmatrix}');

    for (let i = 0; i < Ks.length; i++) {
        let rowString;
        if (i % 2 === 0) {
            rowString = "\\delta_" + (i / 2);
        } else {
            rowString = "\\phi_" + ((i - 1) / 2);
        }
        if (u_zero.includes(i)) {
            rowString = '0';
        }
        if (i < Ks.length) {
            matrixString = matrixString.concat(rowString, '\\\\');
        }
    }

    matrixString = matrixString.concat('\\end{bmatrix} = \\begin{bmatrix}');

    // Reactions
    for (let i = 0; i < Ks.length; i++) {
        let rowString;
        if (i % 2 === 0) {
            rowString = "R_{s," + (i / 2) + '}';
        } else {
            rowString = "M_{s," + ((i - 1) / 2) + '}';
        }
        if (R_zero.includes(i)) {
            rowString = '0';
        }
        if (i < Ks.length) {
            matrixString = matrixString.concat(rowString, '\\\\');
        }
    }

    matrixString = matrixString.concat('\\end{bmatrix} + \\begin{bmatrix}');

    // External forces
    for (let i = 0; i < Ks.length; i++) {
        const F_n_str = F_n.map(x => Math.round(x * 100) / 100);
        const rowString = F_n_str[i].toString() + "\\\\";
        matrixString = matrixString.concat(rowString);
    }

    matrixString = matrixString.concat('\\end{bmatrix}$$');

    p.innerHTML = matrixString;
    MathJax.Hub.Queue(["Typeset", MathJax.Hub, p]);
}
