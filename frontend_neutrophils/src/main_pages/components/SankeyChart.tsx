import * as React from "react"; 
import * as d3 from 'd3';
import * as ReactDOM from "react-dom";

import * as d3sankey from "d3-sankey";

export interface D3SankeyChartProps { id: string, graph: {nodes: any, links: any}, graphtitle?:string }
export interface D3SankeyChartState { }


export default class D3SankeyChart extends React.Component<D3SankeyChartProps, D3SankeyChartState> {

    constructor(props) {
        super(props);
    }
   
    componentWillMount()
    {

    }

    componentWillReceiveProps(newprops)
    {
        this.drawChart(newprops);
    }

    force = null;
    svg = null;

    drawChart(theprops) {

        var self = this;

        var domNode = ReactDOM.findDOMNode(this.refs.graph) as HTMLElement;

        while (domNode.firstChild)
        {
            domNode.removeChild(domNode.firstChild);
        }


        var margin = {top:0, left:0, bottom:0, right:0 };      
        var width = domNode.clientWidth;
        var height = domNode.clientHeight;

        this.svg = d3.select(domNode).append("svg")
        .attr("width", width)
        .attr("height", height);
        //.call(d3.zoom().on("zoom", function () {
        //    self.svg.attr("transform", d3.event.transform)
        //}));
        //.append("g");
        //.attr("transform","translate("+width/2.0+","+height/2.0+")");

        var formatNumber = d3.format(",.0f"),
            format = function (d: any) { return formatNumber(d) + " TWh"; },
            color = d3.scaleOrdinal(d3.schemeCategory10);

        var sankey = d3sankey.sankey()
            .nodeWidth(15)
            .nodePadding(10)
            .extent([[1, 1], [width - 1, height - 6]]);

        sankey(this.testdata);

        var link = this.svg.append("g")
            .attr("class", "links")
            .attr("fill", "none")
            .attr("stroke", "#000")
            .attr("stroke-opacity", 0.2)
            .selectAll("path");

        var node = this.svg.append("g")
            .attr("class", "nodes")
            .attr("font-family", "sans-serif")
            .attr("font-size", 10)
            .selectAll("g");

            link = link
            .data(this.testdata.links)
            .enter().append("path")
            .attr("d", d3sankey.sankeyLinkHorizontal())
            .attr("stroke-width", function (d: any) { return Math.max(1, d.width); });

        link.append("title")
            .text(function (d: any) { return d.source.name + " → " + d.target.name + "\n" + format(d.value); });

        node = node
            .data(this.testdata.nodes)
            .enter().append("g");

        node.append("rect")
            .attr("x", function (d: any) { return d.x0; })
            .attr("y", function (d: any) { return d.y0; })
            .attr("height", function (d: any) { return d.y1 - d.y0; })
            .attr("width", function (d: any) { return d.x1 - d.x0; })
            .attr("fill", function (d: any) { return color(d.name.replace(/ .*/, "")); })
            .attr("stroke", "#000");

        node.append("text")
            .attr("x", function (d: any) { return d.x0 - 6; })
            .attr("y", function (d: any) { return (d.y1 + d.y0) / 2; })
            .attr("dy", "0.35em")
            .attr("text-anchor", "end")
            .text(function (d: any) { return d.name; })
            .filter(function (d: any) { return d.x0 < width / 2; })
            .attr("x", function (d: any) { return d.x1 + 6; })
            .attr("text-anchor", "start");

        node.append("title")
            .text(function (d: any) { return d.name + "\n" + format(d.value); });
    }

    total_width(d)
    {
        // or add d.evidence+d.predicted
        return 20;
    }

    line_shift(d, direction) {
      var rel_x = d.target.x-d.source.x;
      var rel_y = d.target.y-d.source.y;
      var theta = Math.atan2(rel_y, rel_x);
  
      var theta_perpendicular = theta + (Math.PI / 2)*direction;
  
    /*console.log("Theta:" + theta);
      console.log("Theta Per:" + theta_perpendicular);
      console.log("Sin:" + Math.sin(theta_perpendicular));
      console.log("Cos:" + Math.cos(theta_perpendicular));
      console.log("Width:" + width);*/
  
      var width = this.total_width(d)
      var delta_x = width / 4 * Math.cos(theta_perpendicular) // Both lines are pushed 1/4, making it 1/2 away from each other.
      var delta_y = width / 4 * Math.sin(theta_perpendicular) // Both lines are pushed 1/4, making it 1/2 away from each other.
      return [delta_x, delta_y]
    }

    evidence_link_width(d)
    {
        var baseWidth = Math.min(20.0, d.evidence) / 20.0;
        return 5.0 + 5.0*baseWidth;
    }

    prediction_link_width(d)
    {
        var baseWidth = Math.min(20.0, d.predicted) / 20.0;
        return 5.0 + 5.0*baseWidth;
    }

    dragstarted(d) {
        if (!d3.event.active) this.force.alphaTarget(0.5).restart();
        d.fx = d.x;
        d.fy = d.y;
    }

    dragended(d) {
        if (!d3.event.active) this.force.alphaTarget(0.5);
        //d.fx = null;
        //d.fy = null;
    } 

    dragged(d) {
        d.fx = d3.event.x;
        d.fy = d3.event.y;
    }     


    updateGraph()
    {
        console.log("D3 SVG Force parallel updateGraph()");
        console.log(this.props.graph);
    
        if (this.props.graph)
        {
            this.drawChart(this.props);

        }
    }

    willUnmount:boolean = false;

    componentUnmount()
    {
        return this.willUnmount;
    }

    componentWillUnmount()
    {
        this.willUnmount = true;
        console.log("Component Will Unmount");
    }

    componentWillUpdate()
    {
        console.log("Component Will Update");
    }

    componentDidMount()
    {
        console.log("Component Did Mount");
        this.updateGraph();
    }

    componentDidUpdate()
    {
        console.log("D3N4J View DidUpdate")
        this.updateGraph();
    }

    componentDidCatch()
    {
        console.log("Some error occured");
    }

    render() {

        console.log("d3 n4j did render");

        var graphtitle = <span></span>;

        if (this.props.graphtitle)
        {
            graphtitle = <h1>{this.props.graphtitle}</h1>;
        }

        //{width: '800px', height:'400px'}
            return (<div>{graphtitle}<div ref="graph" style={{height: "800px", width: "1200px"}}></div></div>);
    }


    testdata = {
        "links": [
            {
                "source": 15,
                "target": 18,
                "value": 8
            },
            {
                "source": 18,
                "target": 25,
                "value": 4
            },
            {
                "source": 25,
                "target": 0,
                "value": 9
            },
            {
                "source": 18,
                "target": 5,
                "value": 4
            },
            {
                "source": 5,
                "target": 35,
                "value": 25
            },
            {
                "source": 25,
                "target": 35,
                "value": 20
            },
            {
                "source": 25,
                "target": 1,
                "value": 2
            },
            {
                "source": 15,
                "target": 27,
                "value": 47
            },
            {
                "source": 27,
                "target": 5,
                "value": 18
            },
            {
                "source": 27,
                "target": 26,
                "value": 13
            },
            {
                "source": 26,
                "target": 1,
                "value": 2
            },
            {
                "source": 26,
                "target": 35,
                "value": 14
            },
            {
                "source": 27,
                "target": 25,
                "value": 15
            },
            {
                "source": 27,
                "target": 12,
                "value": 8
            },
            {
                "source": 12,
                "target": 35,
                "value": 7
            },
            {
                "source": 12,
                "target": 0,
                "value": 5
            },
            {
                "source": 26,
                "target": 0,
                "value": 4
            },
            {
                "source": 27,
                "target": 14,
                "value": 3
            },
            {
                "source": 14,
                "target": 0,
                "value": 6
            },
            {
                "source": 27,
                "target": 2,
                "value": 3
            },
            {
                "source": 2,
                "target": 0,
                "value": 6
            },
            {
                "source": 5,
                "target": 1,
                "value": 1
            },
            {
                "source": 15,
                "target": 30,
                "value": 8
            },
            {
                "source": 30,
                "target": 26,
                "value": 4
            },
            {
                "source": 30,
                "target": 5,
                "value": 1
            },
            {
                "source": 30,
                "target": 25,
                "value": 3
            },
            {
                "source": 30,
                "target": 12,
                "value": 1
            },
            {
                "source": 30,
                "target": 14,
                "value": 2
            },
            {
                "source": 30,
                "target": 2,
                "value": 2
            },
            {
                "source": 15,
                "target": 8,
                "value": 2
            },
            {
                "source": 8,
                "target": 5,
                "value": 2
            },
            {
                "source": 15,
                "target": 16,
                "value": 4
            },
            {
                "source": 16,
                "target": 25,
                "value": 3
            },
            {
                "source": 16,
                "target": 12,
                "value": 2
            },
            {
                "source": 16,
                "target": 26,
                "value": 2
            },
            {
                "source": 16,
                "target": 14,
                "value": 1
            },
            {
                "source": 16,
                "target": 2,
                "value": 1
            },
            {
                "source": 15,
                "target": 21,
                "value": 1
            },
            {
                "source": 21,
                "target": 25,
                "value": 1
            },
            {
                "source": 15,
                "target": 7,
                "value": 4
            },
            {
                "source": 7,
                "target": 5,
                "value": 1
            },
            {
                "source": 7,
                "target": 25,
                "value": 3
            },
            {
                "source": 7,
                "target": 12,
                "value": 2
            },
            {
                "source": 7,
                "target": 26,
                "value": 1
            },
            {
                "source": 25,
                "target": 6,
                "value": 1
            },
            {
                "source": 12,
                "target": 6,
                "value": 1
            },
            {
                "source": 15,
                "target": 34,
                "value": 1
            },
            {
                "source": 34,
                "target": 5,
                "value": 1
            },
            {
                "source": 5,
                "target": 0,
                "value": 1
            },
            {
                "source": 15,
                "target": 13,
                "value": 1
            },
            {
                "source": 13,
                "target": 25,
                "value": 1
            }
        ],
        "nodes": [
            {
                "name": "location_blood"
            },
            {
                "name": "location_tissue"
            },
            {
                "name": "Gelatinase"
            },
            {
                "name": "lysozyme"
            },
            {
                "name": "HBP"
            },
            {
                "name": "mUnknown"
            },
            {
                "name": "homeostasis"
            },
            {
                "name": "T cell"
            },
            {
                "name": "peritoneal macrophage"
            },
            {
                "name": "tertiary granules"
            },
            {
                "name": "BPI"
            },
            {
                "name": "cathepsin"
            },
            {
                "name": "Chemokine"
            },
            {
                "name": "mucosal type mast cell"
            },
            {
                "name": "Metalloproteinase"
            },
            {
                "name": "Kupffer cell"
            },
            {
                "name": "monocyte"
            },
            {
                "name": "secretory vesicles"
            },
            {
                "name": "hepatocyte"
            },
            {
                "name": "histaminase"
            },
            {
                "name": "inflamation_sterile"
            },
            {
                "name": "endothelial cell"
            },
            {
                "name": "defensin"
            },
            {
                "name": "chemokine receptor"
            },
            {
                "name": "inflamation_infection"
            },
            {
                "name": "Cytokine"
            },
            {
                "name": "Myeloperoxidase"
            },
            {
                "name": "neutrophil"
            },
            {
                "name": "CD35"
            },
            {
                "name": "NET"
            },
            {
                "name": "mature neutrophil"
            },
            {
                "name": "Lactoferrin"
            },
            {
                "name": "homeostasis_excl"
            },
            {
                "name": "secondary granules"
            },
            {
                "name": "splenocyte"
            },
            {
                "name": "cUnknown"
            },
            {
                "name": "granules"
            },
            {
                "name": "neutrophil-derived cathelicidin"
            }
        ]
    };

}