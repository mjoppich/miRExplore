import * as React from "react"; 
import * as d3 from 'd3';
import * as ReactDOM from "react-dom";

import * as d3sankey from "d3-sankey";

export interface D3SankeyChartProps {graph: {nodes: any, links: any}, graphtitle?:string }
export interface D3SankeyChartState { }


export default class D3SankeyChart extends React.Component<D3SankeyChartProps, D3SankeyChartState> {

    constructor(props) {
        super(props);
    }

    public static defaultProps: Partial<D3SankeyChartProps> = {
        graphtitle: "",
        graph: {'nodes': [], 'links': []}
    };
   
    componentWillMount()
    {

    }

    componentWillReceiveProps(newprops)
    {
    }

    force = null;
    svg = null;

    drawChart(theprops) {

        console.log("Sankey Chart Props");
        console.log(theprops);


        var self = this;

        var domNode = ReactDOM.findDOMNode(this.refs.graph) as HTMLElement;

        while (domNode.firstChild)
        {
            domNode.removeChild(domNode.firstChild);
        }

        if (!theprops.graph)
        {
            console.log("Sankey Chart Props Graph Empty - testdata");
            theprops.graph = this.testdata;
        }

        console.log("Sankey Chart Graph");
        console.log(theprops.graph);
        console.log(domNode);

        if (theprops.graph.nodes.length == 0)
        {
            console.log("empty drawChart for empty nodes");
            return;
        }

        var margin = {top:0, left:0, bottom:0, right:0 };      
        var width = domNode.clientWidth;
        var height = domNode.clientHeight;

        console.log("Sankey width")
        console.log(width);
        console.log("Sankey height")
        console.log(height);

        this.svg = d3.select(domNode).append("svg")
        .attr("width", width)
        .attr("height", height);
        //.call(d3.zoom().on("zoom", function () {
        //    self.svg.attr("transform", d3.event.transform)
        //}));
        //.append("g");
        //.attr("transform","translate("+width/2.0+","+height/2.0+")");

        var formatNumber = d3.format(",.0f"),
            format = function (d: any) { return formatNumber(d) + " evidences"; },
            color = d3.scaleOrdinal(d3.schemeCategory10);

        var sankey = d3sankey.sankey()
            .nodeWidth(15)
            .nodePadding(10)
            .extent([[1, 1], [width - 1, height - 6]]);

        sankey(theprops.graph);

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
            .data(theprops.graph.links)
            .enter().append("path")
            .attr("d", d3sankey.sankeyLinkHorizontal())
            .attr("stroke-width", function (d: any) { return Math.max(1, d.width); });

        link.append("title")
            .text(function (d: any) { return d.source.name + " → " + d.target.name + "\n" + format(d.value); });

        node = node
            .data(theprops.graph.nodes)
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
            .text(function (d: any) { var obos = "" ? d.obo == null : d.obo.join(", ") + "\n"; return d.name + "\n" + obos + format(d.value); });

        console.log("Sankey chart drawn.")
    }



    willUnmount:boolean = false;

    componentUnmount()
    {
        return this.willUnmount;
    }

    componentWillUnmount()
    {
        this.willUnmount = true;
        console.log("Sankey Chart Component Will Unmount");
    }

    componentWillUpdate(nextProps, nextState)
    {
        console.log("Sankey Chart Component Will Update");
        console.log(nextProps);
        console.log(nextState);
    }

    componentDidMount()
    {
        console.log("Sankey Chart Component Did Mount");
        this.drawChart(this.props); 
    }

    componentDidUpdate()
    {
        console.log("Sankey Chart Component Did Update");
        this.drawChart(this.props);
    }

    render() {

        console.log("Sankey Chart render");

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