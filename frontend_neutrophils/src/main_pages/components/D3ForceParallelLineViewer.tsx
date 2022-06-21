import * as React from "react"; 
import * as d3 from 'd3';
import * as ReactDOM from "react-dom";

import axios from 'axios';



interface SimNode extends d3.SimulationNodeDatum {
    id: string;
    group: number;
    r: number;
}

interface SimLink extends d3.SimulationLinkDatum<SimNode> {
    value: number;
    d: number;
    s: number;
}

interface Graph {
    nodes: Array<SimNode>;
    links: Array<SimLink>;
}

export interface D3ParallelLinesProps { id: string, graph: {nodes: any, links: any} }
export interface D3ParallelLinesState { }


export default class D3ParallelLinesGraph extends React.Component<D3ParallelLinesProps, D3ParallelLinesState> {

    canvas: any = null;

    canvasSimulation: any = null;
    canvasTransform: any = null;
    context: any = null;
    nodeRadius: number = 20;
    graph: {nodes: any, links: any} = {nodes:[], links:[]};

    constructor(props) {
        super(props);

        console.log(this.props);
    }
   
    componentWillMount()
    {

    }

    drawChart() {

        var self = this;

        var domNode = ReactDOM.findDOMNode(this.refs.graph) as HTMLElement;

        while (domNode.firstChild)
        {
            domNode.removeChild(domNode.firstChild);
        }


        var margin = {top:0, left:0, bottom:0, right:0 };      
        var width = domNode.clientWidth;
        var height = domNode.clientHeight;
        
        if (this.canvasSimulation == null)
        {

            console.log("Creating new canvas");

            this.canvasSimulation = d3.forceSimulation<SimNode, SimLink>();
            this.canvasSimulation
            .force("link", d3.forceLink().id(function(d:any) { return d.index; }).distance(100))
            .force("charge", d3.forceManyBody().strength(function (d, i) {
                var a = i == 0 ? -2000 : -1000;
                return a;
            }).distanceMin(200).distanceMax(1000))
            .force("center", d3.forceCenter(width / 2, height / 2))
            .force("y", d3.forceY(0.01))
            .force("x", d3.forceX(0.01))

            /*
            .force("link", d3.forceLink().id(function(d:any) { return d.index; }).distance(100))
            .force("collide", d3.forceCollide( function(d){return self.nodeRadius + 8 }).iterations(16).strength(-10) )
            .force("charge", d3.forceManyBody().strength(-20).distanceMin(100).distanceMax(150))
            .force("center", d3.forceCenter(width / 2, height / 2))
            //.force("force-x", d3.forceX().strength(0.1))
            //.force("force-y", d3.forceY().strength(0.1));
            */
            console.log("width: " + width);
            console.log("height: " + height);

            this.canvasTransform = d3.zoomIdentity;

        this.canvasSimulation.nodes(this.props.graph.nodes).on("tick", () => {this.canvasRender();});
        var d3simforce: d3.ForceLink<SimNode, SimLink> = this.canvasSimulation.force("link");
        d3simforce.links(this.props.graph.links);

        }

       
        var graphDiv = d3.select( domNode );

        this.canvas = graphDiv.append("canvas");
        this.canvas.attr('width', width);
        this.canvas.attr('height', height);
        this.canvas.style.width = width+"px";
        this.canvas.style.height = height+"px";

        this.context = this.canvas.node().getContext("2d", {preserveDrawingBuffer: true});
        
        this.canvas.call(d3.drag().subject( () => self.getDragSubject() ).on("start", () => self.dragNodeStarted() ).on("drag", () => self.canvasDragged()))
        .call(d3.zoom().scaleExtent([1 / 2, 8]).on("zoom", () => this.canvasZoomed()))
        .call(() => this.canvasRender());

    }
    
        
    canvasZoomed() {
        this.canvasTransform = d3.event.transform;
      this.canvasRender();
    }
    
    getDragSubject() {
      var i,
          x = this.canvasTransform.invertX(d3.event.x),
          y = this.canvasTransform.invertY(d3.event.y),
          dx,
          dy,
          mindist = 100;

    
      for (i = this.props.graph.nodes.length - 1; i >= 0; --i) {
        var point = this.props.graph.nodes[i];
        
        dx = x - point.x;
        dy = y - point.y;
        
        var dist = dx*dx + dy*dy;
        
        if (dist < mindist)
        {
        mindist = dist;
        }
        
        if (dist < this.nodeRadius * this.nodeRadius) {
          point.x = this.canvasTransform.applyX(point.x);
          point.y = this.canvasTransform.applyY(point.y);
          
          return point;
        }
        
      }
      
      //console.log("mindist:" + mindist);
      //console.log("ths: " + circRad * circRad);
    }
    
    canvasDragged() {
      d3.event.subject.fx = this.canvasTransform.invertX(d3.event.x);
      d3.event.subject.fy = this.canvasTransform.invertY(d3.event.y);
      this.canvasRender();
    }
    
    canvasRender() {

        //console.log("canvas render");

        if ((this.context == null) ||(this.context == {}))
        {
            console.log("Updating context");
            this.context = this.canvas.node().getContext('2d');

            return;
        }

        if (this.componentUnmount())
        {
            return;
        }

        if (!('save' in this.context))
        {

            console.log("context=");
            console.log(this.context);
    
        }

        //console.log(this.context);

        this.context.save();
        this.context.clearRect(0, 0, this.canvas.attr('width'), this.canvas.attr('height'));
        
        //this.context.beginPath();
        this.context.translate(this.canvasTransform.x, this.canvasTransform.y);
        this.context.scale(this.canvasTransform.k, this.canvasTransform.k);
        
        //this.context.beginPath();
        //this.context.lineWidth = 20;
        //this.context.strokeStyle = "red";
        this.props.graph.links.forEach( (d) => this.drawLink(d, this.context, -1, d.predicted));
        //this.context.stroke();   

        //this.context.beginPath();
        this.props.graph.links.forEach( (d) => this.drawLink(d, this.context, 1, d.evidence));
        //this.context.stroke();
      
        this.context.beginPath();
        this.props.graph.nodes.forEach((d) => this.drawNode(d, this.context));
        this.context.fill();
      
        this.context.beginPath();
        this.props.graph.nodes.forEach((d) => this.drawNodeLabel(d, this.context));
        this.context.fill();
      
      this.context.restore();
    }
    
    drawLabel( ctx, text, p1, p2, alignment, padding ){
              if (!alignment) alignment = 'center';
              if (!padding) padding = 0;
    
              var dx = p2.x - p1.x;
              var dy = p2.y - p1.y;
              var p, pad;

              if (alignment=='center'){
                p = p1;
                pad = 1/2;
              } else {
                var left = alignment=='left';
                p = left ? p1 : p2;
                pad = padding / Math.sqrt(dx*dx+dy*dy) * (left ? 1 : -1);
              }
              
              var addRotate = 0;
              
              if (dx < 0)
              {
              addRotate += Math.PI
              }
    
              ctx.save();
              ctx.textAlign = alignment;
              ctx.translate(p.x+dx*pad,p.y+dy*pad);
              ctx.rotate(Math.atan2(dy,dx)+addRotate);
              ctx.fillStyle = "black";
              ctx.fillText(text, 0,-2);
              ctx.restore();
            }
    

    
        drawLink(d, ctx, shiftFactor, lineWidth) {
    
            // evidence
            var sourceShift = this.line_shift(d, shiftFactor);

            var startX = d.source.x-sourceShift[0];
            var startY = d.source.y-sourceShift[1];

            var stopX = d.target.x-sourceShift[0];
            var stopY = d.target.y-sourceShift[1];

            ctx.beginPath();

            if (shiftFactor == 1)
            {
                ctx.lineWidth = d.evidence;
                ctx.strokeStyle = "green";
            } else {
                ctx.lineWidth = d.predicted;
                ctx.strokeStyle = "red";
            }

            ctx.moveTo(startX, startY);
            ctx.lineTo(stopX, stopY);

            ctx.stroke();

            
            if (shiftFactor != 1)
            {
                this.drawLabel(ctx, "Hallo", d.source, d.target, null, null);
            }
            
    
          }

          total_width(d)
          {
              // or add d.evidence+d.predicted
              return 40;
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
    
          drawNode(d, ctx) {
    
          
            //ctx.moveTo(d.x - 0.5*this.nodeRadius, d.y-0.5*this.nodeRadius);
            ctx.moveTo(d.x, d.y);
            ctx.fillStyle = "gray";
            ctx.arc(d.x, d.y, this.nodeRadius, 0*Math.PI,2*Math.PI);
            ctx.fill();
            
          }
          
          drawNodeLabel(d, ctx) {
                  ctx.save();
            ctx.textAlign = "center";
            ctx.font = "20px Arial";
            ctx.fillStyle="blue";
            ctx.fillText(d.id, d.x, d.y); 
            ctx.restore();
            }
    
    
          
          
          dragNodeStarted() {
    
            if (!d3.event.active) this.canvasSimulation.alphaTarget(0.3).restart();
    
            
            d3.event.subject.fx = this.canvasTransform.invertX(d3.event.subject.x);
            d3.event.subject.fy = this.canvasTransform.invertY(d3.event.subject.y);
    
          }
    
    
    
    
          dragNodeEnded() {
    
            if (!d3.event.active) this.canvasSimulation.alphaTarget(0);
            d3.event.subject.fx = null;
            d3.event.subject.fy = null;
    
          }


    updateGraph()
    {
        console.log("D3 n4j viewer updateGraph()");

        //var d3graph = this.neo4jgraph2d3(this.props.graph);

        console.log(this.props.graph);
    
        this.drawChart();
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

        //{width: '800px', height:'400px'}
            return (<div><p>graph</p><div ref="graph" style={{height: "500px", width: "1000px"}}></div></div>);
    }
}