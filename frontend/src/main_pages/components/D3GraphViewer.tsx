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

export interface D3Neo4JViewerProps { id: string }
export interface D3Neo4JViewerState { graph: any }


export default class D3Neo4JViewer extends React.Component<D3Neo4JViewerProps, D3Neo4JViewerState> {

    canvas: any = null;

    canvasSimulation: any = null;
    canvasTransform: any = null;
    context: any = null;
    nodeRadius: number = 20;
    graph: {nodes: any, links: any} = {nodes:[], links:[]};

    constructor(props) {
        super(props);

        console.log("d3 n4j viewer constr");
        console.log(this.props);
    }
   
    componentWillMount()
    {
   
        axios.get('http://localhost:3000/api/query')
        .then((response)=>{
          
            console.log("Data received");
            console.log(response.data);
            this.setState({ graph: response.data });

          })
        .catch((err)=>{

            console.log("Fetching data failed");
            console.log(err);

            this.setState({ graph: {nodes: [], links: []} });
        })
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

        this.canvasSimulation.nodes(this.state.graph.nodes).on("tick", () => {this.canvasRender();});
        var d3simforce: d3.ForceLink<SimNode, SimLink> = this.canvasSimulation.force("link");
        d3simforce.links(this.state.graph.links);

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

        console.log("get subject with " + this.state.graph.nodes.length);
    
      for (i = this.state.graph.nodes.length - 1; i >= 0; --i) {
        var point = this.state.graph.nodes[i];
        
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

        if ((this.context == null) ||(this.context == {}))
        {
            console.log("Updatng context");
            this.context = this.canvas.node().getContext('2d');

            return;
        }

        if (!('save' in this.context))
        {

            console.log("context=");
            console.log(this.context);
    
        }

      this.context.save();
      this.context.clearRect(0, 0, this.canvas.attr('width'), this.canvas.attr('height'));
    
      this.context.beginPath();
      this.context.translate(this.canvasTransform.x, this.canvasTransform.y);
      this.context.scale(this.canvasTransform.k, this.canvasTransform.k);
      
      this.state.graph.links.forEach( (d) => this.drawLink(d, this.context));
          this.context.strokeStyle = "#aaa";
          this.context.stroke();
      
      
      
          this.context.beginPath();
          this.state.graph.nodes.forEach((d) => this.drawNode(d, this.context));
      this.context.fill();
      
      
      this.context.beginPath();
      this.state.graph.nodes.forEach((d) => this.drawNodeLabel(d, this.context));
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
              ctx.fillText("EDGE",0,-2);
              ctx.restore();
            }
    
    
        drawLink(d, ctx) {
    
            ctx.moveTo(d.source.x, d.source.y);
            ctx.lineTo(d.target.x, d.target.y);
            
            this.drawLabel(ctx, "Hallo", d.source, d.target, null, null);
    
          }
    
    
    
          drawNode(d, ctx) {
    
          
            ctx.moveTo(d.x + this.nodeRadius, d.y);
            ctx.fillStyle = "black";
            ctx.arc(d.x, d.y, this.nodeRadius, 0, 2 * Math.PI);
            
          }
          
          drawNodeLabel(d, ctx) {
                  ctx.save();
            ctx.textAlign = "center";
            ctx.font = "10px Arial";
            ctx.fillStyle="#aaa";
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

        //var d3graph = this.neo4jgraph2d3(this.state.graph);
    
        this.drawChart();
    }

    componentDidMount()
    {
        if (this.state != null)
        {
            //this.updateGraph();
        }
    }

    componentDidUpdate()
    {
        console.log("D3N4J View DidUpdate")
        if (this.state != null)
        {
            this.updateGraph();
        } 
    }

    componentDidCatch()
    {
        console.log("Some error occured");
    }

    render() {

        console.log("d3 n4j did render");

        //{width: '800px', height:'400px'}
            return (<div><p>graph</p><div ref="graph" style={{height: "500px"}}></div></div>);
    }
}