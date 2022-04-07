//user inputs
var E = 1;
var I = 1;
var nodesX = [0,1,2]
var loads = [0.5,0]
var BC = [2,1,1] //0 = free, 1 = pinned, 2 = fixed

//calculate and initiate some stuff
var numElements = 2+(nodesX.length-1)*2;
//var u_zero = new Array(2+(nodesX.length-1)*2).fill(0).map((x,i) => x+i);
//var R_zero = u_zero;
//calculate length of segments & total length
var L = [];
for (let j = 0;j<nodesX.length-1;j++) {
    L[j] = nodesX[j+1] - nodesX[j]
}
var totLength = L.reduce((current,previous) => current+previous);

var u_zero = [];
var R_zero = [];
//loop over nodes to get BCs
for (let i = 0; i<BC.length;i++) {
    if (BC[i]==0) {
    R_zero.push(i*2) //reaction vertical
    R_zero.push(i*2+1) //reaction moment
    } else if (BC[i]==1) {
    u_zero.push(i*2) //disp vertical
    R_zero.push(i*2+1) //reaction moment
    } else if (BC[i]==2) {
    u_zero.push(i*2) //disp vertical
    u_zero.push(i*2+1) //disp rotation
    }
}
geoDraw(nodesX,loads,BC)


var KeAll = [];
for (let i=0;i<nodesX.length-1;i++) {
KeAll[i] = [
    [12,6*L[i],-12,6].map(x => x/L[i]**3),
    [6*L[i],4*L[i]**2,-6*L[i],2*L[i]**2].map(x => x/L[i]**3),
    [-12,-6*L[i],12*L[i],-6*L[i]].map(x => x/L[i]**3),
    [6*L[i],2*L[i]**2,-6*L[i],4*L[i]**2].map(x => x/L[i]**3)
];
}

function compile(KeAll) {
    var Ks = [];
    sz = 2+KeAll.length*2;
    for (let i=0; i<sz; i++) {
        Ks[i] = new Array(sz).fill(0);
    }

    for (let i=0;i<KeAll.length;i++) {
        id_r = i*2;
        id_c = i*2;
        for (let j=0; j<4;j++) {
        Ks[id_r+j][id_c+0] = Ks[id_r+j][id_c+0] + KeAll[i][j][0]
        Ks[id_r+j][id_c+1] = Ks[id_r+j][id_c+1] + KeAll[i][j][1]
        Ks[id_r+j][id_c+2] = Ks[id_r+j][id_c+2] + KeAll[i][j][2]
        Ks[id_r+j][id_c+3] = Ks[id_r+j][id_c+3] + KeAll[i][j][3]
        }
    }
return Ks
}

var Ks = compile(KeAll)

var u = []
var R = []
for (let i = 0;i<Ks.length;i++) {
    u[i] = 'x' + (i+1)
    R[i] = 'y' + (i+1)
}

var F_n = [];
for (let i = 0; i<nodesX.length;i++) {
    if (i==0) {
        F_n[2*i] = -loads[i]*L[0]/2;
        F_n[2*i+1] = -loads[i]*L[0]**2/12;
    } else if (i==(nodesX.length-1)) {
        F_n[i*2] = -loads[i-1]*L[i-1]/2;
        F_n[i*2+1] = loads[i-1]*L[i-1]**2/12;    
    } else {
        F_n[i*2] = -loads[i-1]*L[i-1]/2 - loads[i]*L[i]/2;
        F_n[i*2+1] = loads[i-1]*L[i-1]**2/12 - loads[i]*L[i]**2/12;
    }
}

//eq1 = Ks[0][0]

//unwrap matrix to system of equations
eq = [];
for (let j = 0; j<sz; j++) {
eq.push(
    Ks[j].reduce((x,y,i) => 
    x.toString() + '+' + y.toString() + '*' + u[i], 
    R[j] + '+' + F_n[j] + '=')
)
}

//add boundary conditions
for (let i = 0; i<u_zero.length; i++) {
    eq.push(u[u_zero[i]] + '= 0')
}
for (let i = 0; i<R_zero.length; i++) {
    eq.push(R[R_zero[i]] + '= 0')
}

//solve system of linear equations with nerdamer library
    var sol2 = nerdamer.solveEquations(eq,
        [...u,...R]);
    console.log(sol2.toString());

function geoDraw(nodesX,loads,BC) {

    var path = [];
    for (let j=0;j<nodesX.length;j++) {
        x = nodesX[j]
        if (BC[j]==1) {
        path[j] = 'M' + x + ',0L' + (x+0.1) + ',-0.1L' + (x-0.1) + ',-0.1 Z' 
        } else if (BC[j]==2) {
        path[j] = 'M' + x + ',0.1L' + x + ',-0.1L' + (x-0.1) + ',-0.1L' + (x-0.1) + ',0.1 Z' 
        } else if (BC[j]==0) {
        path[j] = 'M' + x + ',0L0,0 Z' 
        }
    }

var shapes = [];
for (let j=0;j<nodesX.length;j++) {
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
for (let i=0;i<nodesX.length-1;i++) {
    loadShapes[i] = {
        type:'rect',
        x0:nodesX[i],
        x1:nodesX[i+1],
        y0:0.1,
        y1:loads[i]/2+0.1,
        fillcolor:'rgba(0,0.5,0,0.5)',
        line: {width:1},
    }
}
shapes = [...shapes, ...loadShapes]
var layout = { shapes:shapes,
 xaxis: {scaleanchor: "y",range:[-0.1,0.1+nodesX[nodesX.length-1]]},
 yaxis: {visible:false}
 }
 var trace1 = {
     x:nodesX,
     y:new Array(nodesX.length).fill(0),
     line: {color: 'rgb(0,0,0)',width:5},
     marker: {size:10, line:{color:'rgb(255,0,0)',width:2,fillcolor:'rgb(255,255,255)'}}
    };
var data = [trace1]
//plotly
var myAx = document.getElementById("myFig")
Plotly.newPlot(myAx,data,layout)
}
geoDraw(nodesX,loads,BC)


function linspace(min,max,nel) {
    var nel = nel-1
    var xMax = max; //max x   
    var xMin = min; //min x from previous section
    var x = [...Array(nel+1).keys()];
    x = x.map(a => a * ((xMax - xMin) / nel) + xMin)
    return x
  }

//get reactions
var Rv = [];
for (let i = 0;i<nodesX.length;i++) {
    id = (nodesX.length*2)+i*2
    Rv[i] = sol2[id][1]
}

var xRv = nodesX
var F = []; //points loads from left to right
var xF = [];
var x = linspace(0,totLength,1000);
var y = linspace(0,0,1000);

//first all reactions
for (let i=0;i<(Rv.length);i++) {
    y = x.map((a,j) => a>=xRv[i]?y[j]+Rv[i]:y[j])    
}
//second all UDL
dx = x[2]-x[1];
x_start_udl = [0,1,3] //start of udl
x_end_udl = [1,3,4]
var q = loads
y_udl = linspace(0,0,1000)
for (let i=0;i<x_start_udl.length;i++) {
    for (let j=1;j<x.length;j++) {
        if (x[j]>=x_start_udl[i] && x[j]<=x_end_udl[i]) {
        y_udl[j] = y_udl[j-1] - dx*q[i]
    } else {
        y_udl[j] = Math.min(y_udl[j],y_udl[j-1]);
    }
}
}

y_tot = y.map((y,i) => y+y_udl[i])

var M_tot = []

var M_0 = sol2[sz+1][1];
M_tot[0] = -M_0;
for (let i=1;i<y_tot.length;i++) {
    M_tot[i] = M_tot[i-1]+dx*(y_tot[i]+y_tot[i+1])/2
}

ax = document.getElementById('myPlot');
ax2 = document.getElementById('myPlot2');
trace = {
    x:x,
    y:y_tot,
    mode:'line',
    fill: 'tozeroy'
}
layout = {
    xaxis: {range:[-0.1,totLength+0.1], title:""},
    yaxis: {range:[-1,1], title:"V(x)"},
}
data = [trace]
Plotly.newPlot(ax,data,layout)

var trace_M = {
    x:x,
    y:M_tot.map(x=>x*5),
    mode:'line',
    fill: 'tozeroy'
}
layout2 = {
    xaxis: {range:[-0.1,totLength+0.1], title:""},
    yaxis: {range:[-1,1], title:"M(x)"},
}
var data2 = [trace_M]
Plotly.newPlot(ax2,data2,layout2)