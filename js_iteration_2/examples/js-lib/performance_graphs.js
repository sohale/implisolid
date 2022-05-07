'use strict';

var time_queue = [];

var time_queue__last_timeout = 0;

function report_time(time_msec, callback){
    time_queue.push(time_msec);
    if(time_queue__last_timeout){
        clearTimeout(time_queue__last_timeout);
        time_queue__last_timeout = 0;
    }
    time_queue__last_timeout = setTimeout(function(){
        console.log("Time: (msec)");
        console.log(time_queue);
        if(callback)
            callback();
        time_queue = [];
        time_queue__last_timeout =0;
    }, 1000);
}


var global_histogram_accumulator = [];
function hist_delayed() {
    var a = time_queue;
    //var q = [];
    if(global_histogram_accumulator.length==0){
        global_histogram_accumulator.push(["Time"]);
    }
    var UPPER_BOUND = 12;
    for(var i=0;i<a.length;i++){
        if(a[i]<=UPPER_BOUND)
            global_histogram_accumulator.push([a[i]*1.]); //The precision of the timing is 1 msec :(
    }
    show_hist(global_histogram_accumulator);
}

function hist(){
    //hist_delayed();
}

//setTimeout(function()}{
function add_graph_div(){
    //<div id="chart_div" style="width: 900px; height: 500px;"></div>
    //var elemDiv = document.createElement('div');
    //document.body.appendChild(elemDiv);
    //not nice:
    //document.body.innerHTML += "<div id=\"chart_div\" style=\"width: 900px; height: 500px;\"></div>";
}

'use strict';
var ggraph = 0;
function show_hist(data_array){
if(ggraph==0){
    add_graph_div();
    google.charts.load("current", {packages:["corechart"]});
}
google.charts.setOnLoadCallback(drawChart);
function drawChart() {
    var data = google.visualization.arrayToDataTable(data_array);
    var options = {
        title: 'Time ( Milliseconds )',
        legend: { position: 'none' },
        colors: ['#e7711c'],
        histogram: {
            bucketSize: 1,
            //maxNumBuckets: 200,
            //minValue: -1,
            //maxValue: 200
        },
        //chartArea: { width: 401 },
        hAxis: {
          //ticks: [-1, 0,1,2,3,4,5,6,7,8,9,10,12,14,20,30,40,50,60,80,100,120,140,160,180,200],
          //ticks: [ -1, 0,1,2,3,4,5,6,7,8,9,10,11,12 ],
          //Array.from(Array(10).keys())
          //[...Array(10).keys()]
          ticks: (Array.from(Array(30*1.).keys())).map(function(x){return x/1.;}, Number),
        },
        //vAxis: { scaleType: 'mirrorLog' }
    };
    //see https://developers.google.com/chart/interactive/docs/gallery/histogram#controlling-buckets
    var div_element = document.getElementById('chart_div');
    console.log(div_element);
    var chart = new google.visualization.Histogram(div_element);
    chart.draw(data, options);
}
ggraph++;
}
