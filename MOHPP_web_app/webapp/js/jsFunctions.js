import "./eel.js";
import "https://cdn.plot.ly/plotly-2.4.2.min.js";

async function getVal() {
  let value = await eel.SensorsValues()();
   console.log( value);
   plotval('plot', value);
   request_values();
   
  }
 
           
eel.expose(request_values);
function request_values(){
  getVal();
            	
}

function plotVal(div, value){
  
  const data = []
  var theta = Math.PI/2 
  for (let i = 0; i<4;i++){
      var x = [0,value[i]*Math.cos(Math.PI*2 +theta)];
      var y = [0,value[i]*Math.sin(Math.PI*2 +theta)]; 
      data[data.length] = {x, y, type:'scatter'};
      theta -=Math.PI/2; 
  }


var layout = {
grid: {
    rows: 1,
    columns: 1,
    pattern: 'independent',
    roworder: 'bottom to top'}
};

Plotly.newPlot('myDiv', data, layout);
}

eel.getVal();