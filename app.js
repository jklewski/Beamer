nbeamsSelector = document.querySelector("#beamSelector")
nbeamsSelector.addEventListener("change", function () {
    inputs = document.getElementsByClassName("type1")
    nbeams = parseInt(nbeamsSelector.value);
    for (let i = 0; i < inputs.length; i++) {
        if (i < nbeams) {
            inputs[i].style.display = 'block'
        } else if (i >= nbeams) {
            inputs[i].style.display = 'none'
        }
    }
})

ploadsSelector = document.getElementById("ploadSelector")
ploadsSelector.addEventListener("change", function () {
    inputs = document.getElementsByClassName("type2")
    npoints = parseInt(ploadSelector.value);
    for (let i = 0; i < inputs.length; i++) {
        if (i < npoints) {
            inputs[i].style.display = 'block'
        } else if (i >= npoints) {
            inputs[i].style.display = 'none'
        }
    }
})

$(".lengths, .loads, .ploads,.ploadsx,.custom-select").change(mainFunction)

//run mainfunction on startup
mainFunction()

function mainFunction() {

    //get inputs
    var lengths = document.getElementsByClassName("lengths")
    var inLoads = document.getElementsByClassName("loads")
    var inPLoads = document.getElementsByClassName("ploads")
    var inPLoadsx = document.getElementsByClassName("ploadsx")

    var L = [];
    var loads = [];
    //GET NODES
    for (let i = 0; i < lengths.length; i++) {
        if (lengths[i].parentElement.parentElement.parentElement.style.display.includes('block')) {
            L.push(parseFloat(lengths[i].value))
        }
        if (inLoads[i].parentElement.parentElement.parentElement.style.display.includes('block')) {
            loads.push(parseFloat(inLoads[i].value) / 100)
        }
    }
    var nodesX = [0]
    for (let i = 0; i < L.length; i++) {
        nodesX[i + 1] = nodesX[i] + L[i]
    }

    A = document.getElementsByClassName('ploadsx')
    for (let i = 0; i < A.length; i++) {
        A[i].max = nodesX[nodesX.length - 1]
    }

    var pointLoads = { x: [], y: [] };
    A = []
    for (let i = 0; i < inPLoads.length; i++) {
        if (inPLoads[i].parentElement.parentElement.parentElement.style.display.includes('block')) {
            pointLoads.y.push(parseFloat(inPLoads[i].value / 100))
            pointLoads.x.push(parseFloat(inPLoadsx[i].value))
        }
    }



    //GET P-LOADS

    //var pointLoads = {x:[0.5],y:[0.1]};

    //GET END SUPPORT CONDITIONS
    var BC = new Array(nodesX.length).fill(1)
    BC[0] = parseInt(document.getElementById('LeftSupportCondition').value)
    BC[BC.length - 1] = parseInt(document.getElementById('RightSupportCondition').value)

    //calculate length of segments & total length
    var L = [];
    for (let j = 0; j < nodesX.length - 1; j++) {
        L[j] = nodesX[j + 1] - nodesX[j]
    }
    var totLength = L.reduce((current, previous) => current + previous);

    var u_zero = [];
    var R_zero = [];
    //loop over nodes to get BCs
    for (let i = 0; i < BC.length; i++) {
        if (BC[i] == 0) {
            R_zero.push(i * 2) //reaction vertical
            R_zero.push(i * 2 + 1) //reaction moment
        } else if (BC[i] == 1) {
            u_zero.push(i * 2) //disp vertical
            R_zero.push(i * 2 + 1) //reaction moment
        } else if (BC[i] == 2) {
            u_zero.push(i * 2) //disp vertical
            u_zero.push(i * 2 + 1) //disp rotation
        }
    }
    geoDraw(nodesX, loads, BC)


    var KeAll = [];
    for (let i = 0; i < nodesX.length - 1; i++) {
        KeAll[i] = [
            [12, 6 * L[i], -12, 6 * L[i]].map(x => x / L[i] ** 3),
            [6 * L[i], 4 * L[i] ** 2, -6 * L[i], 2 * L[i] ** 2].map(x => x / L[i] ** 3),
            [-12, -6 * L[i], 12, -6 * L[i]].map(x => x / L[i] ** 3),
            [6 * L[i], 2 * L[i] ** 2, -6 * L[i], 4 * L[i] ** 2].map(x => x / L[i] ** 3)
        ];
    }

    function compile(KeAll) {
        var Ks = [];
        sz = 2 + KeAll.length * 2;
        for (let i = 0; i < sz; i++) {
            Ks[i] = new Array(sz).fill(0);
        }

        for (let i = 0; i < KeAll.length; i++) {
            id_r = i * 2;
            id_c = i * 2;
            for (let j = 0; j < 4; j++) {
                Ks[id_r + j][id_c + 0] = Ks[id_r + j][id_c + 0] + KeAll[i][j][0]
                Ks[id_r + j][id_c + 1] = Ks[id_r + j][id_c + 1] + KeAll[i][j][1]
                Ks[id_r + j][id_c + 2] = Ks[id_r + j][id_c + 2] + KeAll[i][j][2]
                Ks[id_r + j][id_c + 3] = Ks[id_r + j][id_c + 3] + KeAll[i][j][3]
            }
        }
        return Ks
    }

    var Ks = compile(KeAll)

    //round values for perforamnce (small error expected)
    for (let i = 0; i < Ks.length; i++) {
        Ks[i] = Ks[i].map(x => Math.round(x * 10000000) / 10000000)
    }

    var u = []
    var R = []
    for (let i = 0; i < Ks.length; i++) {
        u[i] = 'x' + (i + 1)
        R[i] = 'y' + (i + 1)
    }

    var F_n = [];
    for (let i = 0; i < nodesX.length; i++) {
        if (i == 0) {
            F_n[2 * i] = -loads[i] * L[0] / 2;
            F_n[2 * i + 1] = -loads[i] * L[0] ** 2 / 12;
        } else if (i == (nodesX.length - 1)) {
            F_n[i * 2] = -loads[i - 1] * L[i - 1] / 2;
            F_n[i * 2 + 1] = loads[i - 1] * L[i - 1] ** 2 / 12;
        } else {
            F_n[i * 2] = -loads[i - 1] * L[i - 1] / 2 - loads[i] * L[i] / 2;
            F_n[i * 2 + 1] = loads[i - 1] * L[i - 1] ** 2 / 12 - loads[i] * L[i] ** 2 / 12;
        }
    }
    //add points loads
    for (let j = 0; j < pointLoads.x.length; j++) {
        var id_end = nodesX.map((x, i) => x >= pointLoads.x[j] ? i : null).filter(x => x /= null)[0];
        var id_start = id_end - 1;
        //horizontal
        x_1 = pointLoads.x[j] - nodesX[id_start]
        x_2 = (nodesX[id_end] - nodesX[id_start]) - x_1;
        Li = x_1 + x_2;
        F_n[2 * id_start] -= pointLoads.y[j] * (3 * x_1 + x_2) * x_2 ** 2 / Li ** 3 //reaction left (negative)
        F_n[2 * id_start + 1] -= pointLoads.y[j] * x_2 ** 2 * x_1 / Li ** 2//moment left (negative)
        F_n[2 * id_end] -= pointLoads.y[j] * (3 * x_2 + x_1) * x_1 ** 2 / Li ** 3 //reaction right (negative)
        F_n[2 * id_end + 1] += pointLoads.y[j] * x_1 ** 2 * x_2 / Li ** 2 //moment right (positive)
    }


    //unwrap matrix to system of equations
    eq = [];
    for (let j = 0; j < sz; j++) {
        eq.push(
            Ks[j].reduce((x, y, i) =>
                x.toString() + '+' + y.toString() + '*' + u[i],
                R[j] + '+' + F_n[j] + '=')
        )
    }

    //add boundary conditions
    for (let i = 0; i < u_zero.length; i++) {
        eq.push(u[u_zero[i]] + '= 0')
    }
    for (let i = 0; i < R_zero.length; i++) {
        eq.push(R[R_zero[i]] + '= 0')
    }

    //solve system of linear equations with nerdamer library
    var sol2 = nerdamer.solveEquations(eq,
        [...u, ...R]);
    console.log(sol2.toString());

    function geoDraw(nodesX, loads, BC) {

        var path = [];
        for (let j = 0; j < nodesX.length; j++) {
            x = nodesX[j]
            if (BC[j] == 1) { //pinned (triangle)
                path[j] = 'M' + x + ',0L' + (x + 0.1) + ',-0.1L' + (x - 0.1) + ',-0.1 Z'
            } else if (BC[j] == 2) { //rigid (rectangle)
                if ((j + 1) < nodesX.length) {
                    path[j] = 'M' + x + ',0.1L' + x + ',-0.1L' + (x - 0.1) + ',-0.1L' + (x - 0.1) + ',0.1 Z'
                }
                else {
                    path[j] = 'M' + x + ',0.1L' + x + ',-0.1L' + (x + 0.1) + ',-0.1L' + (x + 0.1) + ',0.1 Z'
                }
            } else if (BC[j] == 0) {
                path[j] = 'M' + x + ',0L0,0 Z'
            }
        }

        var shapes = [];
        for (let j = 0; j < nodesX.length; j++) {
            shapes[j] =
            {
                type: 'path',
                path: path[j],
                fillcolor: 'rgba(200, 200, 200, 0.5)',
                line: {
                    color: 'rgb(0, 0, 0)',
                    width: 0.5,
                }
            }

        }
        //plot UDL loads
        var loadShapes = [];
        for (let i = 0; i < nodesX.length - 1; i++) {
            loadShapes[i] = {
                type: 'rect',
                x0: nodesX[i],
                x1: nodesX[i + 1],
                y0: 0.1,
                y1: loads[i] / 3 + 0.1,
                fillcolor: 'rgba(0,0.5,0,0.4)',
                line: { width: 0 },
            }
        }

        ploadAnnotations = [];
        for (let i = 0; i < pointLoads.x.length; i++) {
            ploadAnnotations[i] = {
                x: pointLoads.x[i],
                y: 0.01,
                xref: 'x',
                yref: 'y',
                text: 'P<sub>' + (i + 1) + '</sub>',
                font: {
                    size: 25,
                    color: '#ff0000'
                },
                arrowcolor: '#ff0000',
                showarrow: true,
                arrowhead: 9,
                ax: 0,
                ay: -50 - (pointLoads.y[i] * 100)
            }
        }

        shapes = [...shapes, ...loadShapes]
        layout = {
            annotations: ploadAnnotations,
            shapes: shapes,
            xaxis: { scaleanchor: "y", range: [-0.1, totLength + 0.1] },
            yaxis: { visible: true },
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            margin: { l: 0, b: 0, t: 10 },
        }
        var trace1 = {
            x: nodesX,
            y: new Array(nodesX.length).fill(0),
            line: { color: 'rgb(0,0,0)', width: 5 },
            marker: { size: 6, line: { color: 'rgb(0,0,0)', width: 2, fillcolor: 'rgb(255,255,255)' } }
        };
        var data = [trace1]
        //plotly
        var myAx = document.getElementById("myFig")
         
        Plotly.newPlot(myAx, data, layout,{responsive: true})
    }
    geoDraw(nodesX, loads, BC)


    function linspace(min, max, nel) {
        var nel = nel - 1
        var xMax = max; //max x   
        var xMin = min; //min x from previous section
        var x = [...Array(nel + 1).keys()];
        x = x.map(a => a * ((xMax - xMin) / nel) + xMin)
        return x
    }



    //Time to calculate the shear force and moment diagram from reactions and loads!
    x_tot = [];
    y_tot = [];
    V_tot = [];
    M_tot = [];
    v_tot = [];
    v_right = [sol2[1][1]];



    //calculate cumulative reaction force, left to right
    var Ri = [sol2[sz][1]];
    for (let i = 1; i < (sz / 2); i++) {
        Ri[i] = Ri[i - 1] + sol2[sz + i * 2][1]; //right side of beam 1
    }


    //calculate node shear force beam-wise from left side by adding up reaction forces and loads
    //then draw V(x) and integrate to get M(x) and v(x) 
    var Vi = []
    for (let i = 0; i < nodesX.length - 1; i++) {

        var M = [];
        x_start = nodesX[i];
        x_end = nodesX[i + 1];
        n_el = 10000;
        x = linspace(x_start, x_end, n_el)
        dx = (x_end - x_start) / n_el


        //find point loads before reaction
        //pID = pointLoads.x.map((x, i) => (x > x_start && x <= x_end) ? i : null) 
        //pID = pID.filter(x => x != null)//IDs
        pID = pointLoads.x.map((x, i) => (x < x_start) ? i : null)
        pID = pID.filter(x => x != null)//IDs
        if (pID.length > 0) {
            pSumLeft = pID.map(x => pointLoads.y[x]).reduce((curr, prev) => curr + prev)
        } else {
            pSumLeft = 0
        }

        qSumLeft = 0;
        for (let k = 0; k < i; k++) {
            qSumLeft += loads[k] * (nodesX[k+1]-nodesX[k]);
        }
        Vi[i] = Ri[i] - qSumLeft - pSumLeft;


        //initiate array at support shear force
        y = new Array(x.length).fill(Vi[i])

        pID = pointLoads.x.map((x, i) => (x >= x_start && x < x_end) ? i : null)
        pID = pID.filter(x => x != null)//IDs
        //loop over point loads in segment and adjust array y
        for (let m = 0; m < pID.length; m++) {
            y = x.map((n, k) => n >= pointLoads.x[pID[m]] ? y[k] - pointLoads.y[pID[m]] : y[k])
        }

        //loop over x to add UDL
        y = y.map((x, u) => x - loads[i] * dx * u)
        //integrate over shear force to get moment
        for (let j = 1; j < x.length; j++) {
            if (i == 0) {
                M[0] = -sol2[sz + 1][1] //moment at left reaction
            } else {
                M[0] = M_tot[M_tot.length - 1] //moment from left side of support
            }
            M[j] = M[j - 1] + dx * y[j]
        }
        //integrate M twice to get deflection
        a = [];
        v = [];

        a = new Array(x.length).fill(null)
        v = new Array(x.length).fill(null)
        a[0] = sol2[i * 2 + 1][1]
        v[0] = sol2[i * 2][1]
        for (let j = 1; j < x.length; j++) {
            a[j] = a[j - 1] + dx * M[j]
            v[j] = v[j - 1] + dx * a[j]
        }

        x_tot.push(...x);
        y_tot.push(...y);
        M_tot.push(...M);
        v_tot.push(...v);
        v_right.push(v[v.length - 1])
    }
    //check if significant error in deflection

    ax2 = document.getElementById('myPlot2');
    //scale factors
    const k_M = Math.max(...[Math.max(...M_tot), Math.abs(Math.min(...M_tot))])
    const k_V = Math.max(...[Math.max(...y_tot), Math.abs(Math.min(...y_tot))])
    const k_v = Math.max(...[Math.max(...v_tot), Math.abs(Math.min(...v_tot))])

    var error = v_right.map((v, i) => Math.abs((v - sol2[i * 2][1])) / k_v).some(x => x > 0.05)

    Math.abs(...M_tot)
    trace = {
        x: x_tot,
        y: y_tot.map(x => x / k_V),
        mode: 'line',
        fill: 'tozeroy',
        name: 'Shear force, scale: ' + k_V.toExponential(2),
    }
    var trace_M = {
        x: x_tot,
        y: M_tot.map(x => x / k_M),
        mode: 'line',
        fill: 'tozeroy',
        name: 'Moment, scale: ' + k_M.toExponential(2),
    }
    var trace_v = {
        x: x_tot,
        y: v_tot.map(x => x / k_v),
        mode: 'line',
        name: 'Deflection, scale: ' + k_v.toExponential(2),
    }

    annotations = [{
        xref: 'paper',
        yref: 'paper',
        x: 0.5,
        xanchor: 'center',
        y: 0.9,
        yanchor: 'bottom',
        text: 'ups... large error detected',
        showarrow: false,
        font: { color: 'red', size: 30 },
    }]
    if (!error || error) {
        annotations[0].opacity = 0
    }

    layout2 = {
        xaxis: { range: [-0.1, totLength + 0.1], title: "" },
        yaxis: { range: [-1, 1], title: "c(x)"},
        showlegend: true,
        legend: {
            x: 0.95,
            y: 0.05,
            xanchor:'right',
            yanchor:'bottom',
            traceorder: 'normal',
            font: {
              family: 'sans-serif',
              size: 12,
              color: '#000'
            },
            bgcolor: '#E2E2E2',
            bordercolor: '#000000',
            borderwidth: 1
          },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        margin: { l: 0, b: 100, t: 0 },
        annotations: annotations,
    }


    var config = {responsive: true} 
    var data2 = [trace_M, trace, trace_v]
    Plotly.newPlot(ax2, data2, layout2,config)

    //write matrix
    var p = document.getElementById('modal-body-text')
    matrixString = '$$EI\\begin{bmatrix}';
    for (let i = 0; i < Ks.length; i++) {
        rowString = Ks[i].map(x => Math.round(x)).toString()
        rowString = rowString.replaceAll(',', '&')
        if (i < Ks.length) {
            matrixString = matrixString.concat(rowString, '\\\\')
        }
    }
    matrixString = matrixString.concat('\\end{bmatrix}', '\\begin{bmatrix}')
    for (let i = 0; i < Ks.length; i++) {
        if (i % 2 == 0) {
            rowString = "\\delta_" + (i / 2)
        }
        if (!(i % 2 == 0)) {
            rowString = "\\phi_" + ((i - 1) / 2)
        }
        if (u_zero.includes(i)) {
            rowString = '0'
        }
        if (i < Ks.length) {
            matrixString = matrixString.concat(rowString, '\\\\')
        }
    }
    matrixString = matrixString.concat('\\end{bmatrix} = \\begin{bmatrix}')

    //Reactions
    for (let i = 0; i < Ks.length; i++) {
        if (i % 2 == 0) {
            rowString = "R_{s," + (i / 2) + '}'
        }
        if (!(i % 2 == 0)) {
            rowString = "M_{s," + ((i - 1) / 2) + '}'
        }
        if (R_zero.includes(i)) {
            rowString = '0'
        }
        if (i < Ks.length) {
            matrixString = matrixString.concat(rowString, '\\\\')
        }
    }
    matrixString = matrixString.concat('\\end{bmatrix} + \\begin{bmatrix}')
    //external forces
    for (let i = 0; i < Ks.length; i++) {
        F_n_str = F_n.map(x => Math.round(x * 100) / 100)
        rowString = F_n_str[i].toString() + "\\\\"
        matrixString = matrixString.concat(rowString)
    }
    matrixString = matrixString.concat('\\end{bmatrix}$$')
    p.innerHTML = matrixString
    MathJax.Hub.Queue(["Typeset", MathJax.Hub, p])



}


