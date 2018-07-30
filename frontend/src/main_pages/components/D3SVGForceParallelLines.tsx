import * as React from "react"; 
import * as d3 from 'd3';
import * as ReactDOM from "react-dom";


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

export interface D3SVGParallelLinesProps { graph: {nodes: any, links: any}, graphInfo: any, graphtitle?:string }
export interface D3SVGParallelLinesState { }


export default class D3SVGParallelLinesGraph extends React.Component<D3SVGParallelLinesProps, D3SVGParallelLinesState> {

    canvas: any = null;

    canvasSimulation: any = null;
    canvasTransform: any = null;
    context: any = null;
    nodeRadius: number = 20;

    constructor(props) {
        super(props);

        console.log(this.props);
    }
   
    componentWillMount()
    {

    }

    componentWillReceiveProps(newprops)
    {
        this.drawChart(newprops);
    }

    svg = null;
    force = null;

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
        .attr("height", height)
        .call(d3.zoom().on("zoom", function () {
            self.svg.attr("transform", d3.event.transform)
        }))
        .append("g")
        .attr("transform","translate("+width/2.0+","+height/2.0+")");

        var attractForce = d3.forceManyBody().strength(100).distanceMax(500).distanceMin(300);
        var repelForce = d3.forceManyBody().strength(-100).distanceMax(250).distanceMin(20);

        var attractForce = d3.forceManyBody().strength(80).distanceMax(400).distanceMin(80);
        var collisionForce = d3.forceCollide(12).strength(1).iterations(100);                       
                            

        this.force = d3.forceSimulation<SimNode, SimLink>() 
            //.force("charge", d3.forceManyBody().strength(-700).distanceMin(200).distanceMax(500)) 
            .alphaDecay(0.03)
            .velocityDecay(0.2)
            .force("link", d3.forceLink().id(function(d:any) { return d.index }).distance(200))
            //.force("x", d3.forceX(width / 2).strength(0.015))
            //.force("y", d3.forceY(height / 2).strength(0.015))
            .force("repel", d3.forceManyBody().strength(-100).distanceMax(250).distanceMin(20))
            .force("collide", d3.forceCollide().radius(function(d)
            {
                return 20;
            }).iterations(5));
/*
            .force("attractForce",attractForce).force("collisionForce",collisionForce)
            .force("link", d3.forceLink().id(function(d:any) { return d.index }).distance(200))
            .force("center", d3.forceCenter(width / 2, height / 2))
            .force("y", d3.forceY(0.001))
            .force("x", d3.forceX(0.001))
            */

        var color = function (group) {
            if (group == 'gene') {
                return "#aaaaaa"
            } else if (group == 'lncrna') {
                return "#fbc280"
            } else {
                return "#405275"
            }
        };
        
        this.force
        .nodes(theprops.graph.nodes) 
        .force("link").links(theprops.graph.links);

        this.maxGroup1Value = 0;
        this.maxGroup2Value = 0;

        for (var i=0; i < theprops.graph.links.length; ++i)
        {
            var edge = theprops.graph.links[i];

            if (edge.group1 > this.maxGroup1Value)
            {
                this.maxGroup1Value = edge.group1;
            }

            if (edge.group2 > this.maxGroup2Value)
            {
                this.maxGroup2Value = edge.group2;
            }
        }


    var group1Link = this.svg.selectAll(".g1-link")
        .data(theprops.graph.links)
        .enter()
        .append("line")
        .attr("class", "link");

    var group2Link = this.svg.selectAll(".g2-link")
        .data(theprops.graph.links)
        .enter()
        .append("line")
        .attr("class", "link");

    var node = this.svg.selectAll(".node")
        .data(theprops.graph.nodes)
        .enter().append("g")
        .attr("class", "node")
        .call(d3.drag()
        .on("start", this.dragstarted.bind(this))
        .on("drag", this.dragged.bind(this))
        .on("end", this.dragended.bind(this))); 

        node.append('circle')
        .attr('r', 13)
        .attr('fill', function (d) {
            return color(d.group);
        });


    node.append("text")
        .attr("dx", 20)
        .attr("dy", ".35em")
        .style("font-size", "18px")
        .text(function (d) {
            return d.name
        });


        var tip;
        self.svg.on("click", function(){
          if (tip) tip.remove();
        });
        node.on("click", function(d){
          d3.event.stopPropagation(); 
        
          if (tip) tip.remove();
          
          tip  = self.svg.append("g")
            .attr("transform", "translate(" + d.x  + "," + d.y + ")");
            
          var rect = tip.append("rect")
            .style("fill", "white")
            .style("stroke", "steelblue");
          
          tip.append("text")
            .text("Name: " + d.name)
            .attr("dy", "1em")
            .attr("x", 5);
            
          tip.append("text")
            .text("Info: " + d.info)
            .attr("dy", "2em")
            .attr("x", 5);
      
          var con = self.props.graph.links
            .filter(function(d1){
              return d1.source.id === d.id;
            })
            .map(function(d1){
              return d1.target.name + " with weight " + d1.weight;
            })
            
          tip.append("text")
            .text("Connected to: " + con.join(","))
            .attr("dy", "3em")
            .attr("x", 5);
          
          var bbox = tip.node().getBBox();
          rect.attr("width", bbox.width + 5)
              .attr("height", bbox.height + 5)
        });


    group1Link
    .style("stroke-width", function stroke(d)  {return self.group1_link_width(d) })
    .style("stroke", "#70C05A")
    //.style("stroke-width", "3.5px")
    .style("stroke-opacity", "1.0")
  
   group2Link
    .style("stroke-width", function stroke(d)  {return self.group2_link_width(d) })
    .style("stroke", "#438DCA")
    //.style("stroke-width", "3.5px")
    .style("stroke-opacity", "1.0")


    var tickFunction = function () {

        group1Link
        .attr("x1", function(d) { return d.source.x-self.line_shift(d,1)[0]; })
        .attr("y1", function(d) { return d.source.y-self.line_shift(d,1)[1]; })
        .attr("x2", function(d) { return d.target.x-self.line_shift(d,1)[0]; })
        .attr("y2", function(d) { return d.target.y-self.line_shift(d,1)[1]; });
        group2Link
            .attr("x1", function(d) { return d.source.x-self.line_shift(d,-1)[0]; })
            .attr("y1", function(d) { return d.source.y-self.line_shift(d,-1)[1]; })
            .attr("x2", function(d) { return d.target.x-self.line_shift(d,-1)[0]; })
            .attr("y2", function(d) { return d.target.y-self.line_shift(d,-1)[1]; });

        node.attr("transform", function (d) {
            return "translate(" + d.x + "," + d.y + ")";
        });
    };

    this.force.on("tick", tickFunction);

    var     padding = 1.5, // separation between same-color circles
    clusterPadding = 6, // separation between different-color circles
    maxRadius = 12;

    // Resolves collisions between d and all other circles.
    function collide(alpha) {
        var quadtree = d3.quadtree(node);
        return function(d) {
        var r = d.radius + maxRadius + Math.max(padding, clusterPadding),
            nx1 = d.x - r,
            nx2 = d.x + r,
            ny1 = d.y - r,
            ny2 = d.y + r;
        quadtree.visit(function(quad: any, x1, y1, x2, y2) {
            if (quad.point && (quad.point !== d)) {
            var x = d.x - quad.point.x,
                y = d.y - quad.point.y,
                l = Math.sqrt(x * x + y * y),
                r = d.radius + quad.point.radius + (d.cluster === quad.point.cluster ? padding : clusterPadding);
            if (l < r) {
                l = (l - r) / l * alpha;
                d.x -= x *= l;
                d.y -= y *= l;
                quad.point.x += x;
                quad.point.y += y;
            }
            }
            return x1 > nx2 || x2 < nx1 || y1 > ny2 || y2 < ny1;
        });
        };
    }


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

    maxGroup1Value: number;
    maxGroup2Value: number;

    group1_link_width(d)
    {
        var baseWidth = Math.max(20.0, d.group1) / this.maxGroup1Value;
        return 5.0 + 5.0*baseWidth;
    }

    group2_link_width(d)
    {
        var baseWidth = Math.max(20.0, d.group2) / this.maxGroup2Value;
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

    oldDraw()
    {


        /**
         * 
         * 
         * 
         * 
         * OLD
         * 
         * 
         */
        var width = 0;
        var height = 0;
        
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

        var self=this;
        var domNode = null;
       
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
        this.props.graph.links.forEach( (d) => this.drawLink(d, this.context, -1, d.group1));
        //this.context.stroke();   

        //this.context.beginPath();
        this.props.graph.links.forEach( (d) => this.drawLink(d, this.context, 1, d.group2));
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
                ctx.lineWidth = d.group1;
                ctx.strokeStyle = "green";
            } else {
                ctx.lineWidth = d.group2;
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
}