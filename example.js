// with much help from https://observablehq.com/@deckerlukas/line-graph-with-moving-average
// and the function below https://observablehq.com/@d3/moving-average

// compute moving average
function movingAverage(values, N) {
    let i = 0;
    let sum = 0;
    const means = new Float64Array(values.length).fill(NaN);
    for (let n = Math.min(N - 1, values.length); i < n; ++i) {
        sum += values[i];
    }
    for (let n = values.length; i < n; ++i) {
        sum += values[i];
        means[i] = sum / N;
        sum -= values[i - N + 1];
    }
    return means;
}

// set up margins and height

const margin = ({ top: 20, bottom: 120, left: 40, right: 20 })
const margin1 = ({ top: 300, bottom: 40, left: 40, right: 20 })

const height1 = 300
const height2 = 80

const width = 954

const svg = d3
    .select("#vis")
    .append("svg")
    .attr("viewBox", [0, 0, 1000, 400])
    .attr("width", 900)
    .attr("height", 400);

// parse the CSV
const data = d3.csv("./Athaliana_genome/Athaliana_genome_stats.csv", function (d) {

    return {
        ID: d.ID,
        bin: +d.bin.replace(/ [0-9]+-/, ""),
        GCPercent: +d['GC%'],
        GCSkew: +d.GCSkew,
        UniqueKmers: +d.UniqueKmers
    }


})

data.then(function (data) {
    // what are the unique fasta headers?
    const chromosomes = [...new Set(data.map(d => d.ID))]
    const variables = Object.keys(data[0]).slice(2, 5) // the keys of the data

    // populate the chromosome dropdown
    d3.select("#chromosomeDropdown")
        .selectAll('myOptionsChrom')
        .data(chromosomes)
        .join('option')
        .text(d => d) // text showed in the menu
        .attr("value", d => d) // corresponding value returned by the button

    // populate the variable dropdown
    d3.select("#VariableDropdown")
        .selectAll('myOptionsVar')
        .data(variables)
        .join('option')
        .text(d => d) // text showed in the menu
        .attr("value", d => d) // corresponding value returned by the button

    // do a little wrangling
    // test on the first chromosome
    let filteredData = d3.group(data, d => d.ID).get(chromosomes[0])
    let filteredDataMA = movingAverage(filteredData.map(d => d[variables[0]]), 100) // work out moving average input later
    for (let i = 0; i < filteredData.length; i++) {
        filteredData[i].MA = filteredDataMA[i];
    }

    // scales and axes
    let xMin = 0;
    let xMax = d3.max(filteredData.map(d => d.bin));
    let yMin = 0;
    let yMax = d3.max(filteredData.map(d => d[variables[0]]))

    const x = d3
        .scaleLinear()
        .domain([xMin, xMax])
        .range([0, width])

    const y = d3
        .scaleLinear()
        .domain([yMin, yMax])
        .range([height1, 0])

    const x2 = d3
        .scaleLinear()
        .domain([xMin, xMax])
        .range([0, width])

    const y2 = d3
        .scaleLinear()
        .domain([yMin, yMax])
        .range([height2, 50])

    const xAxis = d3.axisBottom(x).tickFormat(x => `${x / 1000000}MB`)
    const yAxis = d3.axisLeft(y)
    const xAxis2 = d3.axisBottom(x2).tickFormat(x => `${x / 1000000}MB`)
    const yAxis2 = d3.axisLeft(y2)

    // line functions

    const line = d3
        .line()
        .x(d => x(d.bin))
        .y(d => y(d[variables[0]]))

    const line2 = d3
        .line()
        .x(d => x2(d.bin))
        .y(d => y2(d[variables[0]]))

    const movAvgLine = d3
        .line()
        .defined(d => !isNaN(d.MA))
        .x(d => x(d.bin))
        .y(d => y(d.MA))

    // the plot itself

    const clip = d3
        .selectAll(svg)
        .append('defs')
        .append('svg:clipPath')
        .attr('id', 'clip')
        .append('svg:rect')
        .attr('width', width)
        .attr('height', height1)
        .attr('x', 0)
        .attr('y', 0);

    const linechart = d3
        .selectAll(svg)
        .append('g')
        .attr('class', 'focus')
        .attr('transform', `translate(${margin.left} ${margin.top})`)
        .attr('clip-path', 'url(#clip)');

    const focus = d3
        .selectAll(svg)
        .append('g')
        .attr('class', 'focus')
        .attr('transform', `translate(${margin.left} ${margin.top})`);

    const context = d3
        .selectAll(svg)
        .append('g')
        .attr('class', 'context')
        .attr('transform', `translate(${margin.left} ${margin1.top})`);

    const focusXaxis = focus
        .append('g')
        .attr('class', 'axis axis--x')
        .attr('transform', `translate(0, ${height1})`)
        .call(xAxis);

    focus
        .append('g')
        .attr('class', 'axis axis--y')
        .call(yAxis);

    // y axis label
    d3.selectAll(svg)
        .append('text')
        .attr("font-family", "sans-serif")
        .attr("font-size", 10)
        .attr('x', 30)
        .attr('y', 13)
        .text('GC%');

    // x axis label
    d3.selectAll(svg)
        .append('text')
        .attr("font-family", "sans-serif")
        .attr("font-size", 10)
        .attr('x', width - margin.right - 35)
        .attr('y', height1 + 13)
        .text(data.x);

    const mainLine = linechart
        .append('path')
        .datum(filteredData)
        .attr('class', 'line')
        .attr('d', line)
        .style('fill', 'none')
        .attr('stroke', "black")
        .attr('stroke-width', 1);

    const movingAverageLine = linechart
        .append('path')
        .datum(filteredData)
        .attr('class', 'line1')
        .attr('d', movAvgLine)
        .style('fill', 'none')
        .attr('stroke', 'red')
        .attr('stroke-width', 0.8);

    const bottomLine = context
        .append('path')
        .datum(filteredData)
        .attr('class', 'line2')
        .attr('d', line2)
        .style('fill', 'none')
        .attr('stroke', "black")
        .attr('stroke-width', 0.1);

    context
        .append('g')
        .attr('class', 'axis axis--x')
        .attr('transform', `translate(0, ${height2})`)
        .call(xAxis2);

    const brush = d3
        .brushX()
        .extent([[0, 40], [width, height2]])
        .on('brush end', function (event, d) {
            if (event.sourceEvent && event.sourceEvent.type === "zoom") return; // ignore brush-by-zoom
            var s = event.selection || x2.range();
            x.domain(s.map(x2.invert, x2));
            linechart.selectAll(".line").attr("d", line);
            linechart.selectAll('.line1').attr('d', movAvgLine);
            linechart.selectAll(".line2").attr("d", line2);
            focus.selectAll(".axis--x").call(xAxis);
        });

    context
        .append('g')
        .attr('class', "brush")
        .call(brush);


    // add here the variable of interest
    function updateChart(movingAverageBin, chromosome, variable) {

        // update the data
        let filteredData = d3.group(data, d => d.ID).get(chromosome)
        let filteredDataMA = movingAverage(filteredData.map(d => d[variable]), movingAverageBin) // work out moving average input later
        for (let i = 0; i < filteredData.length; i++) {
            filteredData[i].MA = filteredDataMA[i];
        }
        // update the x scales
        x.domain([xMin, d3.max(filteredData.map(d => d.bin))])
        x2.domain([xMin, d3.max(filteredData.map(d => d.bin))])
        // update the y scales
        y.domain([variable === "GCSkew" ? d3.min(filteredData.map(d => d[variable])) : 0, d3.max(filteredData.map(d => d[variable]))])
        y2.domain([variable === "GCSkew" ? d3.min(filteredData.map(d => d[variable])) : 0, d3.max(filteredData.map(d => d[variable]))])

        // apply changes to the x axis
        svg.selectAll('.axis--x')
            .transition()
            .duration(1000)
            .call(xAxis)

        // apply the changes to the y axis

        // update the lines
        line.y(d => y(d[variable]))
        line2.y(d => y2(d[variable]))


        // apply changes to the lines
        movingAverageLine
            .datum(filteredData)
            .attr('d', movAvgLine)
        mainLine
            .datum(filteredData)
            .attr('d', line)
        bottomLine
            .datum(filteredData)
            .attr('d', line2)

    }

    // Listen to the slider
    d3.select("#movingAverage").on("change", function (d) {
        selectedValue = this.value
        updateChart(selectedValue, 
            d3.select("#chromosomeDropdown").node().value, 
            d3.select("#VariableDropdown").node().value)
    })
    // Listen to the chromosome dropdown
    d3.select("#chromosomeDropdown").on("change", function (d) {
        selectedGroup = this.value
        updateChart(d3.select("#movingAverage").node().value, 
        selectedGroup, 
        d3.select("#VariableDropdown").node().value)
    })
    // Listen to the variable dropdown
    d3.select("#VariableDropdown").on("change", function (d) {
        selectedGroup = this.value
        updateChart(d3.select("#movingAverage").node().value, 
        d3.select("#chromosomeDropdown").node().value, 
        selectedGroup)
    })

    return svg.node();

})