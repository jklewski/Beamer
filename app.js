//sidebar
const open_btn = document.querySelector(".open-btn")
const close_btn = document.querySelector(".close-btn")
const popup = document.querySelector(".popup")
const main_popup = document.querySelector(".main-popup")

open_btn.addEventListener("click", () => {
    popup.style.display = 'flex';
    main_popup.style.cssText = 'animation: slide-in .5s ease; animation-fill-mode: forwards'
}
)

close_btn.addEventListener('click', () => {
    main_popup.style.cssText = 'animation: slide-out .5s ease; animation-fill-mode: forwards'
    setTimeout(() => {
        popup.style.display = 'none';
    }, 500);
})

window.addEventListener('click', (e) => {
    if (e.target == document.querySelector('.popup-overlay')) {
        main_popup.style.cssText = 'animation: slide-out .5s ease; animation-fill-mode: forwards'
        setTimeout(() => {
            popup.style.display = 'none';
        }, 500);
    }
})


//add some listeners
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

$(".lengths, .loads, .custom-select").change(mainFunction)

//run mainfunction on startup
mainFunction()

function mainFunction() {
    //get inputs
    var lengths = document.getElementsByClassName("lengths")
    var inLoads = document.getElementsByClassName("loads")

    var L = [];
    var loads = [];
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
                fillcolor: 'rgba(255, 140, 184, 0.5)',
                line: {
                    color: 'rgb(255, 140, 184)'
                }
            }

        }
        //plot loads
        var loadShapes = [];
        for (let i = 0; i < nodesX.length - 1; i++) {
            loadShapes[i] = {
                type: 'rect',
                x0: nodesX[i],
                x1: nodesX[i + 1],
                y0: 0.1,
                y1: loads[i] / 2 + 0.1,
                fillcolor: 'rgba(0,0.5,0,0.5)',
                line: { width: 1 },
            }
        }
        shapes = [...shapes, ...loadShapes]
        var layout = {
            shapes: shapes,
            xaxis: { scaleanchor: "y", range: [-0.1, 0.1 + nodesX[nodesX.length - 1]] },
            yaxis: { visible: true, range: [-0.1, 1] },
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            margin: { b: 0, t: 10 },
        }
        var trace1 = {
            x: nodesX,
            y: new Array(nodesX.length).fill(0),
            line: { color: 'rgb(0,0,0)', width: 5 },
            marker: { size: 10, line: { color: 'rgb(255,0,0)', width: 2, fillcolor: 'rgb(255,255,255)' } }
        };
        var data = [trace1]
        //plotly
        var myAx = document.getElementById("myFig")
        Plotly.newPlot(myAx, data, layout)


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

    //get reactions
    var Rv = [];
    for (let i = 0; i < nodesX.length; i++) {
        id = (nodesX.length * 2) + i * 2
        Rv[i] = sol2[id][1]
    }

    var xRv = nodesX
    var x = linspace(0, totLength, 1000);
    var y = linspace(0, 0, 1000);

    //first all reactions
    for (let i = 0; i < (Rv.length); i++) {
        y = x.map((a, j) => a >= xRv[i] ? y[j] + Rv[i] : y[j])
    }

    //second all UDL
    dx = x[2] - x[1];

    x_start_udl = nodesX.slice(0, nodesX.length - 1) //start of udl
    x_end_udl = nodesX.slice(1, nodesX.length)
    var q = loads
    y_udl = linspace(0, 0, 1000)
    for (let i = 0; i < x_start_udl.length; i++) {
        for (let j = 1; j < x.length; j++) {
            if (x[j] >= x_start_udl[i] && x[j] <= x_end_udl[i]) {
                y_udl[j] = y_udl[j - 1] - dx * q[i]
            } else {
                y_udl[j] = Math.min(y_udl[j], y_udl[j - 1]);
            }
        }
    }

    y_tot = y.map((y, i) => y + y_udl[i])

    var M_tot = []

    var M_0 = sol2[sz + 1][1];
    M_tot[0] = -M_0;
    for (let i = 1; i < y_tot.length; i++) {
        M_tot[i] = M_tot[i - 1] + dx * (y_tot[i] + y_tot[i + 1]) / 2
    }

    ax2 = document.getElementById('myPlot2');
    trace = {
        x: x,
        y: y_tot,
        mode: 'line',
        fill: 'tozeroy',
        name: 'Shear force',
    }

    var trace_M = {
        x: x,
        y: M_tot.map(x => x * 5),
        mode: 'line',
        fill: 'tozeroy',
        name: 'Moment',
    }
    layout = {
        xaxis: { range: [-0.1, totLength + 0.1], title: "" },
        yaxis: { range: [-1, 1], title: "c(x)" },
        showlegend: true,
        legend: {
            x: 1,
            xanchor: 'right',
            y: 1
        },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        margin: { b: 100, t: 0 },
    }
    var data2 = [trace_M, trace]
    Plotly.newPlot(ax2, data2, layout)

    //write matrix
    var p = document.getElementById('KMatrix')
    matrixString = '$$\\frac{EI}{L^3}\\begin{bmatrix}';
    for (let i = 0; i < Ks.length; i++) {
        rowString = Ks[i].map(x => Math.round(x)).toString()
        rowString = rowString.replaceAll(',', '&')
        matrixString = matrixString.concat(rowString, '\\\\')
    }
    matrixString = matrixString.concat('\\end{bmatrix}$$')
    p.innerHTML = matrixString
    MathJax.Hub.Queue(["Typeset", MathJax.Hub])

    window.onresize = function () {
        ax = document.getElementsByClassName("responsive-plot")
        for (let i = 0; i < ax.length; i++) {
            ax[i].style.width = "100%"
            Plotly.Plots.resize(ax[i]);
        }

    }

}