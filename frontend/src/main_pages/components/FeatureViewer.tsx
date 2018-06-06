import * as d3 from 'd3';


class ToolTip
{

    tooltipDiv;
    selectedRect;
    bodyNode = null;
    tooltipColor = null;

    obj: any;

    constructor(obj: any, public superThis: FeatureViewer)
    {
        this.bodyNode = this.superThis.d3div.node();
        this.tooltipColor = this.superThis.options.tooltipColor ? this.superThis.options.tooltipColor : "orangered";

        this.obj = obj;

        
    }

    tooltip(selection: any)
    {

        var self = this;

        selection.on('mouseover.tooltip', function (pD, pI) {

            console.log("tooltip mouseover")
            // Clean up lost tooltips
            d3.select('body').selectAll('div.tooltip').remove();
            // Append tooltip
            var absoluteMousePos = d3.mouse(self.bodyNode);
            var rightside = (absoluteMousePos[0] > self.superThis.width);
            if (rightside) {
                self.tooltipDiv = self.superThis.d3div.append('div')
                    .attr('class', 'tooltip3');
            } else {
                self.tooltipDiv = self.superThis.d3div.append('div')
                    .attr('class', 'tooltip2');
                self.tooltipDiv.style("left", absoluteMousePos[0] - 15 + 'px');
           }
            self.tooltipDiv
                .style("bottom", (self.bodyNode.offsetHeight - absoluteMousePos[1] + 16) + 'px')
                .style("background-color", "#eeeeee")
                .style("width", "auto")
                .style("max-width", "170px")
                .style("height", "auto")
                .style("max-height", "68px")
                .style("padding", "5px")
                .style("font", "10px sans-serif")
                .style("text-align", "center")
                .style("position", "absolute")
                .style("z-index", "45")
                .style("box-shadow", "0 1px 2px 0 #656565")
            /*.style({
                bottom: (self.bodyNode.offsetHeight - absoluteMousePos[1] + 16) + 'px',
                'background-color': '#eee',
                width: 'auto',
                'max-width': '170px',
                height: 'auto',
                'max-height': '68px',
                padding: '5px',
                "font": '10px sans-serif',
                'text-align': 'center',
                position: 'absolute',
                'z-index': 45,
                'box-shadow': '0 1px 2px 0 #656565' 
            });*/
            if (self.obj.type === "path") {
                var first_line = '<p style="margin:2px;font-weight:700;color:' + self.tooltipColor +'">' + pD[0].x + '&#x256d;&#x256e;' + pD[1].x + '</p>';
                if (pD.description) var second_line = '<p style="margin:2px;color:' + self.tooltipColor +';font-size:9px">' + pD.description + '</p>';
                else var second_line = '';
            } else if (self.obj.type === "line") {
                var elemHover:any = self.superThis.updateLineTooltip(absoluteMousePos[0],pD);

                if (elemHover.description) {
                    var first_line = '<p style="margin:2px;font-weight:700;color:' + self.tooltipColor +'">' + elemHover.x + ' : <span> ' + elemHover.y + '</span></p>';
                    var second_line = '<p style="margin:2px;color:' + self.tooltipColor +';font-size:9px">' + elemHover.description + '</p>';
                }
                else {
                    var first_line = '<p style="margin:2px;color:' + self.tooltipColor +'">position : <span id="tLineX">' + elemHover.x + '</span></p>';
                    var second_line = '<p style="margin:2px;color:' + self.tooltipColor +'">count : <span id="tLineC">' + elemHover.y + '</span></p>';
                }
            } else if (self.obj.type === "unique" || pD.x === pD.y) {
                var first_line = '<p style="margin:2px;font-weight:700;color:' + self.tooltipColor +'">' + pD.x + '</p>';
                if (pD.description) var second_line = '<p style="margin:2px;color:' + self.tooltipColor +';font-size:9px">' + pD.description + '</p>';
                else var second_line = '';
            } else {
                var first_line = '<p style="margin:2px;font-weight:700;color:' + self.tooltipColor +'">' + pD.x + ' - ' + pD.y + '</p>';
                if (pD.description) var second_line = '<p style="margin:2px;color:' + self.tooltipColor +';font-size:9px">' + pD.description + '</p>';
                else var second_line = '';
            }

            self.tooltipDiv.html(first_line + second_line);
            if (rightside) {
                self.tooltipDiv.style("left", (absoluteMousePos[0] + 10 - (self.tooltipDiv.node().getBoundingClientRect().width)) + 'px');
            }
        })
            .on('mousemove.tooltip', function (pD, pI) {
            
                if (self.obj.type === "line") {
                    var absoluteMousePos = d3.mouse(self.bodyNode);
                    var elemHover: any = self.superThis.updateLineTooltip(absoluteMousePos[0],pD);
                    if (elemHover.description) {
                        var first_line = '<p style="margin:2px;color:' + self.tooltipColor +'">' + elemHover.x + ' : <span> ' + elemHover.y + '</span></p>';
                        var second_line = '<p style="margin:2px;color:' + self.tooltipColor +';font-size:9px">' + elemHover.description + '</p>';
                    }
                    else {
                        var first_line = '<p style="margin:2px;color:' + self.tooltipColor +'">position : <span id="tLineX">' + elemHover.x + '</span></p>';
                        var second_line = '<p style="margin:2px;color:' + self.tooltipColor +'">count : <span id="tLineC">' + elemHover.y + '</span></p>';
                    }
                    self.tooltipDiv.html(first_line + second_line);
//                            $('#tLineX').text(elemHover.x);
//                            $('#tLineC').text(elemHover.y);  
                }
                // Move tooltip
                // IE 11 sometimes fires mousemove before mouseover
                if (self.tooltipDiv === undefined) { return; }
                var absoluteMousePos = d3.mouse(self.bodyNode);
                var rightside = (absoluteMousePos[0] > self.superThis.width);
                if (rightside) {
                    self.tooltipDiv.attr("class", "tooltip3");
                    self.tooltipDiv.style("left", (absoluteMousePos[0] + 10 - (self.tooltipDiv.node().getBoundingClientRect().width)) + 'px')
                                    .style("bottom",(self.bodyNode.offsetHeight - absoluteMousePos[1] + 16) + 'px');
                } else {
                    self.tooltipDiv.attr("class", "tooltip2");
                    self.tooltipDiv.style("left", (absoluteMousePos[0] - 15) + 'px')
                        .style("bottom", (self.bodyNode.offsetHeight - absoluteMousePos[1] + 16) + 'px');
                }
            })
            .on('mouseout.tooltip', function (pD, pI) {
                // Remove tooltip
                self.tooltipDiv.remove();
            })
            .on('click', function (pD, pI) {
                var xTemp;
                var yTemp;
                var xRect;
                var widthRect;
                var elemHover;

                if(this.nodeName === "text") {
                    var rect = document.getElementById(this.previousSibling.id);
                    if(rect != null) self.colorSelectedFeat(rect, self.obj);
                }
                else self.colorSelectedFeat(this, self.obj);

                d3.select('body').selectAll('div.selectedRect').remove();
                // Append tooltip
                self.selectedRect = self.superThis.d3div.append('div')
                    .attr('class', 'selectedRect');
                if (self.obj.type === "path") {
                    xTemp = pD[0].x;
                    yTemp = pD[1].x;
                } else if (self.obj.type === "line") {
                    var absoluteMousePos = d3.mouse(self.bodyNode);
                    elemHover = self.superThis.updateLineTooltip(absoluteMousePos[0],pD);
                    xTemp = elemHover.x - 0.5;
                    yTemp = elemHover.x + 0.5;
                } else if (self.obj.type === "unique" || pD.x === pD.y) {
                    xTemp = pD.x - 0.4;
                    yTemp = pD.y + 0.4;
                } else {
                    xTemp = pD.x;
                    yTemp = pD.y;
                }

                console.log("svgwidth");
                console.log(self.superThis.SVGOptions.brushActive);
                console.log( d3.select(".brush .overlay").attr("width"));
                console.log( self.superThis.svgContainer.node().getBBox().width);

                var svgWidth = self.superThis.SVGOptions.brushActive ? d3.select(".brush .overlay").attr("width") : self.superThis.svgContainer.node().getBBox().width;

                console.log("elem mose click");
                console.log("determining width");
                console.log(xTemp);
                console.log(yTemp);
                console.log(svgWidth);

                if (self.superThis.scaling(xTemp) < 0 && self.superThis.scaling(yTemp) > svgWidth) {
                    xRect = self.superThis.margin.left;
                    widthRect = parseInt(svgWidth) + 5;

                    console.log("determining width1");
                    console.log(widthRect);


                } else if (self.superThis.scaling(xTemp) < 0) {
                    xRect = self.superThis.margin.left;
                    widthRect = (self.superThis.scaling(yTemp));

                    console.log("determining width2");
                    console.log(widthRect);

                } else if (self.superThis.scaling(yTemp) > svgWidth) {
                    xRect = self.superThis.scaling(xTemp) + self.superThis.margin.left;
                    widthRect = parseInt(svgWidth) - self.superThis.scaling(xTemp);
                    widthRect =  widthRect + 5;

                    console.log("determining width3");
                    console.log(widthRect);

                } else {
                    xRect = self.superThis.scaling(xTemp) + self.superThis.margin.left;
                    widthRect = (self.superThis.scaling(yTemp) - self.superThis.scaling(xTemp));

                    console.log("determining width4");
                    console.log(widthRect);
                }
                self.selectedRect
                    .style("left", xRect + "px")
                    .style("top", ($(self.superThis.viewerDiv).find(".svgHeader").length) ? 60 + 'px' : 10 + 'px')
                    .style('background-color', 'rgba(0, 0, 0, 0.2)')
                    .style("width", widthRect + "px")
                    .style("height", (self.superThis.Yposition + 50) + 'px')
                    .style("position", "absolute")
                    .style("z-index", "-1")
                    .style("box-shadow", "0 1px 2px 0 #656565");

                /*
                self.selectedRect.style({
                    left: xRect + 'px',
                    top: ($(self.superThis.div + " .svgHeader").length) ? 60 + 'px' : 10 + 'px',
                    'background-color': 'rgba(0, 0, 0, 0.2)',
                    width: widthRect + 'px',
                    height: (self.superThis.Yposition + 50) + 'px',
                    position: 'absolute',
                    'z-index': -1,
                    'box-shadow': '0 1px 2px 0 #656565'
                });*/

                if (CustomEvent) {
                    var event = new CustomEvent(self.superThis.events.FEATURE_SELECTED_EVENT, {
                        detail: {
                            start: self.obj.type === "path" ? pD[0].x : self.obj.type === "line" ? elemHover.x : pD.x,
                            end: self.obj.type === "path" ? pD[1].x : self.obj.type === "line" ? elemHover.y : pD.y,
                            id: self.obj.type === "path" ? pD[0].id : self.obj.type === "line" ? elemHover.id : pD.id,
                            description:self.obj.type === "path" ? pD[0].description : self.obj.type === "line" ? elemHover.description : pD.description
                        }
                    });
                    self.superThis.svgElement.dispatchEvent(event);
                } else {
                    console.warn("CustomEvent is not defined....");
                }
                if (self.superThis.trigger) self.superThis.trigger(self.superThis.events.FEATURE_SELECTED_EVENT, {
                    start: self.obj.type === "path" ? pD[0].x : self.obj.type === "line" ? elemHover.x : pD.x,
                    end: self.obj.type === "path" ? pD[1].x : self.obj.type === "line" ? elemHover.y : pD.y,
                    id: self.obj.type === "path" ? pD[0].id : self.obj.type === "line" ? elemHover.id : pD.id,
                    description:self.obj.type === "path" ? pD[0].description : self.obj.type === "line" ? elemHover.description : pD.description
                });

            });
    }

    /*
        attr(...args: any[]) {
            if (!arguments.length) return args;
            attrs = _x;
            return this;
        }

        style(_x) {
            if (!arguments.length) return styles;
            styles = _x;
            return this;
        }
        */

        colorSelectedFeat(feat, object) {
            //change color && memorize
            if (this.superThis.featureSelected !== {}) d3.select(this.superThis.featureSelected.id).style("fill", this.superThis.featureSelected.originalColor);
            if (object.type !== "path" && object.type !== "line"){
                this.superThis.featureSelected = {"id": feat, "originalColor": d3.select(feat).style("fill") || object.color};
                d3.select(feat).style("fill", "orangered");
            }
        }
}

class PreComputing
{

    constructor(public superThis: FeatureViewer)
    {

    }

    path(object) {
        object.data.sort(function (a, b) {
            return a.x - b.x;
        });

        var level = this.superThis.addLevel(object.data);


        object.data = object.data.map(function (d) {
            return [{
                x: d.x,
                y: 0,
                id: d.id,
                description: d.description,
                color: d.color
            }, {
                x: d.y,
                y: d.level + 1,
                id: d.id
            }, {
                x: d.y,
                y: 0,
                id: d.id
            }]
        })

        this.superThis.pathLevel = level * 10 + 5;
        object.height = level * 10 + 5;
    }

    line(object) {
        if (!object.height) object.height = 10;
        var shift = parseInt(object.height);
        var level = 0;
        for (var i in object.data) {
            object.data[i].sort(function (a, b) {
                return a.x - b.x;
            });
            if (object.data[i][0].y !== 0) {
                object.data[i].unshift({
                    x:object.data[i][0].x-1,
                    y:0
                })
            }
            if (object.data[i][object.data[i].length -1].y !== 0){
                object.data[i].push({
                    x:object.data[i][object.data[i].length -1].x+1,
                    y:0
                })
            }
            var maxValue = Math.max.apply(Math,object.data[i].map(function(o){return Math.abs(o.y);}));
            level = maxValue > level ? maxValue : level;
            

            object.data[i] = [object.data[i].map(function (d) {
                return {
                    x: d.x,
                    y: d.y,
                    id: d.id,
                    description: d.description
                }
            })]
        }
        this.superThis.lineYscale.range([0, -(shift)]).domain([0, -(level)]);
        this.superThis.pathLevel = shift * 10 +5;
        object.level = level;
        object.shift = shift * 10 +5;
    }

    multipleRect(object) {
        object.data.sort(function (a, b) {
            return a.x - b.x;
        });
        this.superThis.level = this.superThis.addLevel(object.data);
        this.superThis.pathLevel = this.superThis.level * 10 + 5;
    }
};


class FillSVG
{
    preComputing: PreComputing=null;

    constructor(public superThis: FeatureViewer)
    {
        this.preComputing = new PreComputing(superThis);
    }

    typeIdentifier(object: any) {

        var self=this;

        if (object.type === "rect") {
            self.preComputing.multipleRect(object);
            this.superThis.yData.push({
                title: object.name,
                y: this.superThis.Yposition,
                filter: object.filter
            });
            self.rectangle(object, this.superThis.Yposition);
        } else if (object.type === "text") {
            self.sequence(object.data, this.superThis.Yposition);
            self.superThis.yData.push({
                title: object.name,
                y: self.superThis.Yposition,
                filter: object.filter
            });
            self.superThis.scaling.range([5, self.superThis.width-5]);
        } else if (object.type === "unique") {
            self.unique(object, self.superThis.Yposition);
            self.superThis.yData.push({
                title: object.name,
                y: self.superThis.Yposition,
                filter: object.filter
            });
        } else if (object.type === "multipleRect") {
            self.preComputing.multipleRect(object);
            self.multipleRect(object, self.superThis.Yposition, self.superThis.level);
            self.superThis.yData.push({
                title: object.name,
                y: self.superThis.Yposition,
                filter: object.filter
            });
            self.superThis.Yposition += (self.superThis.level - 1) * 10;
        } else if (object.type === "path") {
            self.preComputing.path(object);
            self.path(object, self.superThis.Yposition);
            self.superThis.Yposition += self.superThis.pathLevel;
            self.superThis.yData.push({
                title: object.name,
                y: self.superThis.Yposition - 10,
                filter: object.filter
            });
        } else if (object.type === "line") {
            if (!(Array.isArray(object.data[0]))) object.data = [object.data];
            if (!(Array.isArray(object.color))) object.color = [object.color];
            var negativeNumbers = false;
            object.data.forEach(function(d){
                if (d.filter(function(l){ return l.y < 0}).length) negativeNumbers = true;
            });
            self.preComputing.line(object);
            self.line(object, self.superThis.Yposition);
            self.superThis.Yposition += self.superThis.pathLevel;
            self.superThis.yData.push({
                title: object.name,
                y: self.superThis.Yposition - 10,
                filter: object.filter
            });
            self.superThis.Yposition += negativeNumbers ? self.superThis.pathLevel-5 : 0;
        }
    }

    sequence(seq, position, start=0) {
        //Create group of sequence

        var self=this;

        this.superThis.svgContainer.append("g")
            .attr("class", "seqGroup")
            .selectAll(".AA")
            .data(seq)
            .enter()
            .append("text")
            .attr("clip-path", "url(#clip)")
            .attr("class", "AA")
            .attr("text-anchor", "middle")
            .attr("x", function (d, i) {
                return self.superThis.scaling.range([5, self.superThis.width-5])(i + start)
            })
            .attr("y", position)
            .attr("font-size", "10px")
            .attr("font-family", "monospace")
            .text(function (d, i) {
                return d
            });
    }

    sequenceLine() {
        var self=this;
        //Create line to represent the sequence
        if (self.superThis.SVGOptions.dottedSequence){
            var dottedSeqLine = self.superThis.svgContainer.selectAll(".sequenceLine")
                .data([[{x:1,y:12},{x:self.superThis.fvLength,y:12}]])
                .enter()
                .append("path")
                .attr("clip-path", "url(#clip)")
                .attr("d", self.superThis.line)
                .attr("class","sequenceLine")
                .style("z-index", "0")
                .style("stroke", "black")
                .style("stroke-dasharray","1,3")
                .style("stroke-width", "1px")
                .style("stroke-opacity",0);

            dottedSeqLine
                .transition()
                .duration(500)
                .style("stroke-opacity", 1);
        }
    }

    rectangle(object, position) {

        var self=this;

        //var rectShift = 20;
        if (!object.height) object.height = 12;
        var rectHeight = object.height;
        
        var rectShift = rectHeight + rectHeight/3;
        var lineShift = rectHeight/2 - 6;
//                var lineShift = rectHeight/2 - 6;

        var rectsPro = this.superThis.svgContainer.append("g")
            .attr("class", "rectangle")
            .attr("clip-path", "url(#clip)")
            .attr("transform", "translate(0," + position + ")");
        
        var dataline=[];
        for (var i = 0; i < this.superThis.level; i++) {
            dataline.push([{
                    x: 1,
                    y: (i * rectShift + lineShift)
                }, {
                    x: this.superThis.fvLength,
                    y: (i * rectShift + lineShift)
                }]);
        }
        rectsPro.selectAll(".line" + object.className)
            .data(dataline)
            .enter()
            .append("path")
            .attr("d", self.superThis.line)
            .attr("class", function () {
                return "line" + object.className
            })
            .style("z-index", "0")
            .style("stroke", object.color)
            .style("stroke-width", "1px");


        var rectsProGroup = rectsPro.selectAll("." + object.className + "Group")
            .data(object.data)
            .enter()
            .append("g")
            .attr("class", object.className + "Group")
            .attr("transform", function (d) {
                return "translate(" + self.superThis.rectX(d) + ",0)"
            });

        var objToolTip = new ToolTip(object, self.superThis);

        rectsProGroup
            .append("rect")
            .attr("class", "element " + object.className)
            .attr("id", function (d) {
                return "f" + d.id
            })
            .attr("y", function (d) {
                return d.level * rectShift
            })
            .attr("width", (d) => this.superThis.rectWidth2(d))
            .attr("height", rectHeight)
            .style("fill", function(d) { return d.color || object.color })
            .style("z-index", "13")
            .call((s) => objToolTip.tooltip(s));
            //.call(d3.helper.tooltip(object));

        rectsProGroup
            .append("text")
            .attr("class", "element " + object.className + "Text")
            .attr("y", function (d) {
                return d.level * rectShift + rectHeight/2
            })
            .attr("dy", "0.35em")
            .style("font-size", "10px")
            .text(function (d) {
                return d.description
            })
            .style("fill", "black")
            .style("z-index", "15")
            .style("visibility", (d) => {
                if (d.description) {
                    return (self.superThis.scaling(d.y) - self.superThis.scaling(d.x)) > d.description.length * 8 && rectHeight > 11 ? "visible" : "hidden";
                } else return "hidden";
            })
            .call((s) => objToolTip.tooltip(s));
            //.call(d3.helper.tooltip(object));


        //rectsPro.selectAll("." + object.className)
        //    .data(object.data)
        //    .enter()
        //    .append("rect")
        //    .attr("clip-path", "url(#clip)")
        //    .attr("class", "element "+object.className)
        //    .attr("id", function(d) { return "f"+d.id })
        //    .attr("x", X)
        //    .attr("width", rectWidth)
        //    .attr("height", 12)
        //    .style("fill", object.color)
        //    .style("z-index", "13")
        //    .call(d3.helper.tooltip(object));

        this.superThis.forcePropagation(rectsProGroup);
        var uniqueShift = rectHeight > 12 ? rectHeight - 6 : 0;
        this.superThis.Yposition += this.superThis.level < 2 ? uniqueShift : (this.superThis.level-1) * rectShift + uniqueShift;
    }

    unique(object, position) {
        var rectsPro = this.superThis.svgContainer.append("g")
            .attr("class", "uniquePosition")
            .attr("transform", "translate(0," + position + ")");

        var self=this;

        var dataline=[];
        dataline.push([{
                x: 1,
                y: 0
            }, {
                x: self.superThis.fvLength,
                y: 0
            }]);
        
        rectsPro.selectAll(".line" + object.className)
            .data(dataline)
            .enter()
            .append("path")
            .attr("clip-path", "url(#clip)")
            .attr("d", self.superThis.line)
            .attr("class", "line" + object.className)
            .style("z-index", "0")
            .style("stroke", object.color)
            .style("stroke-width", "1px");

        var objToolTip = new ToolTip(object, self.superThis);
        
        rectsPro.selectAll("." + object.className)
            .data(object.data)
            .enter()
            .append("rect")
            .attr("clip-path", "url(#clip)")
            .attr("class", "element " + object.className)
            .attr("id", function (d) {
                return "f" + d.id
            })
            .attr("x", function (d) {
                return self.superThis.scaling(d.x - 0.4)
            })
            .attr("width", function (d) {
                if (self.superThis.scaling(d.x + 0.4) - self.superThis.scaling(d.x - 0.4) < 2) return 2;
                else return self.superThis.scaling(d.x + 0.4) - self.superThis.scaling(d.x - 0.4);
            })
            .attr("height", 12)
            .style("fill", function(d) {return d.color ||  object.color})
            .style("z-index", "3")
            .call((s) => objToolTip.tooltip(s));

            this.superThis.forcePropagation(rectsPro);
    }

    path(object, position) {

        var self=this;
        var pathsDB = this.superThis.svgContainer.append("g")
            .attr("class", "pathing")
            .attr("transform", "translate(0," + position + ")");

        var dataline=[];
        dataline.push([{
                x: 1,
                y: 0
            }, {
                x: self.superThis.fvLength,
                y: 0
            }]);
        
        pathsDB.selectAll(".line" + object.className)
            .data(dataline)
            .enter()
            .append("path")
            .attr("clip-path", "url(#clip)")
            .attr("d", self.superThis.lineBond)
            .attr("class", "line" + object.className)
            .style("z-index", "0")
            .style("stroke", object.color)
            .style("stroke-width", "1px");
            
        var objToolTip = new ToolTip(object, self.superThis);

        pathsDB.selectAll("." + object.className)
            .data(object.data)
            .enter()
            .append("path")
            .attr("clip-path", "url(#clip)")
            .attr("class", "element " + object.className)
            .attr("id", function (d) {
                return "f" + d[0].id
            })
            .attr("d", self.superThis.lineBond)
            .style("fill", "none")
            .style("stroke", function(d) {return d[0].color || object.color})
            .style("z-index", "3")
            .style("stroke-width", "2px")
            .call((s) => objToolTip.tooltip(s));

            this.superThis.forcePropagation(pathsDB);
    }

    line(object, position) {
        if (!object.interpolation) object.interpolation = "monotone";
        if (object.fill === undefined) object.fill = true;
        var histog = this.superThis.svgContainer.append("g")
            .attr("class", "lining")
            .attr("transform", "translate(0," + position + ")");

        var self=this;

        var dataline=[];
        dataline.push([{
                x: 1,
                y: 0
            }, {
                x: self.superThis.fvLength,
                y: 0
            }]);
        
        histog.selectAll(".line" + object.className)
            .data(dataline)
            .enter()
            .append("path")
            .attr("clip-path", "url(#clip)")
            .attr("d", self.superThis.lineBond)
            .attr("class", "line" + object.className)
            .style("z-index", "0")
            .style("stroke", "black")
            .style("stroke-width", "1px");

            var objToolTip = new ToolTip(object, self.superThis);


        object.data.forEach(function(dd,i,array){
            histog.selectAll("." + object.className + i)
            .data(dd)
            .enter()
            .append("path")
            .attr("clip-path", "url(#clip)")
            .attr("class", "element " + object.className + " " + object.className + i)
            .attr("d", self.superThis.lineGen.interpolate(object.interpolation))
            .style("fill", object.fill ? self.superThis.shadeBlendConvert(0.6, object.color[i]) || self.superThis.shadeBlendConvert(0.6, "#000") : "none")
            .style("stroke", object.color[i] || "#000")
            .style("z-index", "3")
            .style("stroke-width", "2px")
//                    .style("shape-rendering", "crispEdges")
            .call((s) => objToolTip.tooltip(s));
})
        
        this.superThis.forcePropagation(histog);
    }

    multipleRect(object, position, level) {
        var rectHeight = 8;
        var rectShift = 10;
        var rects = this.superThis.svgContainer.append("g")
            .attr("class", "multipleRects")
            .attr("transform", "translate(0," + position + ")");

        for (var i = 0; i < level; i++) {
            rects.append("path")
                .attr("d", this.superThis.line([{
                    x: 1,
                    y: (i * rectShift - 2)
                }, {
                    x: this.superThis.fvLength,
                    y: (i * rectShift - 2)
                }]))
                .attr("class", function () {
                    return "line" + object.className
                })
                .style("z-index", "0")
                .style("stroke", object.color)
                .style("stroke-width", "1px");
        }

        var self=this;

        rects.selectAll("." + object.className)
            .data(object.data)
            .enter()
            .append("rect")
            .attr("clip-path", "url(#clip)")
            .attr("class", "element " + object.className)
            .attr("id", function (d) {
                return "f" + d.id
            })
            .attr("x", self.superThis.X)
            .attr("y", function (d) {
                return d.level * rectShift
            })
            .attr("width", this.superThis.rectWidth)
            .attr("height", rectHeight)
            .style("fill", function(d) { return d.color || object.color })
            .style("z-index", "13")
            .call(new ToolTip(object, self.superThis));

        this.superThis.forcePropagation(rects);
    }
};

class Transition
{
    constructor(public superThis: FeatureViewer)
    {

    }


        rectangle(object) {

            var self=this;

            this.superThis.svgContainer.selectAll(".line" + object.className)
                .attr("d",self.superThis.line.x((d) => this.superThis.scaling(d.x)));
            var transit1, transit2;

            if (this.superThis.animation) {
                transit1 = this.superThis.svgContainer.selectAll("." + object.className + "Group")
                    .transition()
                    .duration(500);
                transit2 = this.superThis.svgContainer.selectAll("." + object.className)
                    .transition()
                    .duration(500);
            }
            else {
                transit1 = this.superThis.svgContainer.selectAll("." + object.className + "Group");
                transit2 = this.superThis.svgContainer.selectAll("." + object.className);
            }
            transit1.attr("transform", function (d) {
                        return "translate(" + self.superThis.rectX(d) + ",0)"
                    });

            transit2
                .attr("width", (d) => this.superThis.rectWidth2(d));
                this.superThis.svgContainer.selectAll("." + object.className + "Text")
                .style("visibility", (d) => {
                    if (d.description) {
                        return (this.superThis.scaling(d.y) - this.superThis.scaling(d.x)) > d.description.length * 8 && object.height > 11 ? "visible" : "hidden";
                    } else return "hidden";
                });
        }

        multiRec(object) {
            var self=this;
            this.superThis.svgContainer.selectAll("." + object.className)
//                    .data(object.data)
                //.transition()
                //.duration(500)
                .attr("x", (d) => {
                    return this.superThis.scaling(d.x)
                })
                .attr("width", (d) => {
                    return this.superThis.scaling(d.y) - this.superThis.scaling(d.x)
                });
        }

        unique(object) {
            var self=this;

            self.superThis.svgContainer.selectAll(".line" + object.className)
                .attr("d",self.superThis.line.x(function (d) {
                return self.superThis.scaling(d.x);
            }));
            var transit;
            if (self.superThis.animation) {
                transit = self.superThis.svgContainer.selectAll("." + object.className)
//                    .data(object.data)
                    .transition()
                    .duration(500);
            }
            else {
                transit = self.superThis.svgContainer.selectAll("." + object.className);
            }
            
            transit
//                    .data(object.data)
                //.transition()
                //.duration(500)
                .attr("x", function (d) {
                    return self.superThis.scaling(d.x - 0.4)
                })
                .attr("width", function (d) {
                    if (self.superThis.scaling(d.x + 0.4) - self.superThis.scaling(d.x - 0.4) < 2) return 2;
                    else return self.superThis.scaling(d.x + 0.4) - self.superThis.scaling(d.x - 0.4);
                });
        }

        path(object) {

            var self=this;

            self.superThis.svgContainer.selectAll(".line" + object.className)
                .attr("d",self.superThis.lineBond.x(function (d) {
                                        return self.superThis.scaling(d.x);
                                    })
                                  .y(function (d) {
                                        return -d.y * 10 + object.height;
                                    })
                     );
            var transit;
            if (self.superThis.animation) {
                transit = self.superThis.svgContainer.selectAll("." + object.className)
//                    .data(object.data)
                    .transition()
                    .duration(500);
            }
            else {
                transit = self.superThis.svgContainer.selectAll("." + object.className);
            }
            transit
                .attr("d", self.superThis.lineBond.y(function (d) {
                    return -d.y * 10 + object.height;
                }));
        }

        line(object) {

            var self=this;

            self.superThis.lineYscale.range([0, -(object.height)]).domain([0, -(object.level)]);
            self.superThis.svgContainer.selectAll(".line" + object.className)
                .attr("d", self.superThis.lineGen.y(function (d) {
                    return self.superThis.lineYscale(-d.y) * 10 + object.shift;
                }));
            var transit;
            if (self.superThis.animation) {
                transit = self.superThis.svgContainer.selectAll("." + object.className)
//                    .data(object.data)
                    .transition()
                    .duration(500);
            }
            else {
                transit = self.superThis.svgContainer.selectAll("." + object.className);
            }
            
            transit
                .attr("d", self.superThis.lineGen.y(function (d) {
                    return self.superThis.lineYscale(-d.y) * 10 + object.shift;
                })
                      .interpolate(object.interpolation)
                     );
        }

        text(object, start) {

            
            var transit;
            if (this.superThis.animation) {
                transit = this.superThis.svgContainer.selectAll("." + object.className)
//                    .data(object.data)
                    .transition()
                    .duration(500);
            }
            else {
                transit = this.superThis.svgContainer.selectAll("." + object.className);
            }
            transit
                .attr("x", (d, i) => this.superThis.scaling(i + start));
            }
}

class MouseEventData implements MouseEventInit {
    public button: number = 0;
    public buttons: number = 0;
    public clientX: number = 0;
    public clientY: number = 0;
    public relatedTarget: EventTarget = null;
    public screenX: number = 0;
    public screenY: number = 0;
    public pageX: number = 0;
    public pageY: number = 0;
}

export default class FeatureViewer {

    transitionObj = new Transition(this);

    events: any;
    //div = null;
    //el = null;
    svgElement;
    sequence = null;
    intLength = 0;
    fvLength = 0;
    features = [];
    SVGOptions = {
        showSequence: false,
        brushActive: false,
        verticalLine: false,
        dottedSequence: true
    };

    isInt(value:string) {

        return Number(value) && Number(value) % 1 === 0;        
    }

    offset = null;

    pathLevel = 0;
    svg;
    svgContainer;
    filter;
    yData = [];
    yAxisSVG;
    yAxisSVGgroup;
    Yposition = 20;
    level = 0;
    seqShift = 0;
    enableZoom = false;
    zoomMax = 50;
    current_extend = null;
    featureSelected = null;
    animation = null;

    scaling = null;



    updateLineTooltip(mouse,pD){
        var xP = mouse-110;
        var elemHover = {};
        for (var l=0; l<pD.length;l++) {
            if (this.scaling(pD[l].x) < xP && this.scaling(pD[l+1].x) > xP) {
                if ((xP - this.scaling(pD[l].x)) < (this.scaling(pD[l+1].x) - xP )) {
                    elemHover = pD[l];
                }
                else elemHover = pD[l+1];
                break;
            }
        }
        return elemHover;
    }

    scalingPosition = null;

    transition_data(features, start) {
        var self = this;

        features.forEach(function (o) {
            if (o.type === "rect") {
                self.transitionObj.rectangle(o);
            } else if (o.type === "multipleRect") {
                self.transitionObj.multiRec(o);
            } else if (o.type === "unique") {
                self.transitionObj.unique(o);
            } else if (o.type === "path") {
                self.transitionObj.path(o);
            } else if (o.type === "line") {
                self.transitionObj.line(o);
            } else if (o.type === "text") {
                self.transitionObj.text(o, start);
            }
        });
    }

    width: number;
    height: number;
    margin: any;
    options: any;

    fillSVG: FillSVG = null;
    tooltip: ToolTip = null;

    viewerDiv: Element = null;
    d3div = null;


    constructor(sequence, div, options) {
//        var nxSeq = sequence.startsWith('NX_') ? true : false;
        var self = this;

        this.fillSVG = new FillSVG(this);


        this.viewerDiv = div;
        this.d3div = d3.select(div);
        //this.el = document.getElementById(div.substring(1));
        //this.el = div;

        this.sequence = sequence;
        this.intLength =  this.isInt(sequence) ? sequence : null;
        this.fvLength = this.intLength | sequence.length;

        // if (!div) var div = window;
        this.events = {
          FEATURE_SELECTED_EVENT: "feature-viewer-position-selected",
            FEATURE_DESELECTED_EVENT: "feature-viewer-position-deselected",
          ZOOM_EVENT: "feature-viewer-zoom-altered"
        };

        // if (!div) var div = window;

        this.offset = {start:1,end:this.fvLength};
        if (options && options.offset) {
            this.offset = options.offset;
            if (this.offset.start < 1) {
                this.offset.start = 1;
                console.warn("WARNING ! offset.start should be > 0. Thus, it has been reset to 1.");
            }
        }

        this.options = options;

        this.current_extend = { 
                    length : this.offset.end - this.offset.start,
                    start : this.offset.start,
                    end : this.offset.end
                }
        this.featureSelected = {};
        this.animation = true;
        /**
         * Private Methods
         */

            //Init box & scaling
        d3.select(div)
            .style("position", "relative")
            .style("padding", "0px")
            .style("z-index", "2");

        this.margin = {
                top: 10,
                right: 20,
                bottom: 20,
                left: 110
            };
            
        this.width = $(div).width() - this.margin.left - this.margin.right - 17,
        this.height = 600 - this.margin.top - this.margin.bottom;

        console.log(this.width);
        console.log(this.height);
        console.log(this.offset);

        this.scaling = d3.scaleLinear()
            .domain([this.offset.start, this.offset.end])
            .range([5, self.width-5]);

        this.scalingPosition = d3.scaleLinear()
            .domain([0, self.width])
            .range([this.offset.start, this.offset.end]);
        
    this.lineYscale = d3.scaleLinear()
        .domain([0,-30])
        .range([0,-20]);

        
    this.line = d3.line()
        .curve(d3.curveLinear)
        .x(function (d) {
            return self.scaling(d[0]);
        })
        .y(function (d) {
            return d[1] + 6;
        });

        this.xAxis = d3.axisBottom( self.scaling ).tickFormat(d3.format("d"));

    //Create Axis
    /*
    this.xAxis = d3.svg.axis()
        .scale(self.scaling)
        .tickFormat(d3.format("d"))
        .orient("bottom");
*/

        this.yAxisScale = d3.scaleBand()
        .domain(["0", ""+self.yData.length])
        .rangeRound([0, 500])
        .padding(0.1);

        this.yAxis = d3.axisLeft(self.yAxisScale).tickValues(self.yData).tickFormat(function (d) {
            return d.toString();
        });



            /*
        this.yAxis = d3.svg.axis()
            .scale(self.yAxisScale)
            .tickValues(self.yData) //specify an array here for values
            .tickFormat(function (d) {
                return d
            })
            .orient("left");
*/

            this.brush = d3.brushX()
            .extent([[self.scaling.range()[0], 0], [self.scaling.range()[1], self.height]])
            //.x(self.scaling)
            //.on("brush", brushmove)
            .on("end", () => this.brushend());            


            this.lineBond = d3.line()
            .curve(d3.curveStepBefore)
            .x(function (d) {
                return self.scaling(d[0]);
            })
            .y(function (d) {
                return -d[1] * 10 + self.pathLevel;
            });
    
            this.lineGen = d3.line()
    //          .interpolate("cardinal")
          .x(function(d) {
            return self.scaling(d[0]);
          })
          .y(function (d) {
                return self.lineYscale(-d[1]) * 10 + self.pathLevel;
            });

            $(window).on("resize", this.resizeCallback.bind(this));

            this.initSVG(this.viewerDiv, options);

        }

        lineBond: any;
        brush: any;
        xAxis: any;
        yAxis: any;
        yAxisScale: any;
        lineGen: any;
        lineYscale: any;
        line: any;

        //COMPUTING FUNCTION
        X(d) {
            return this.scaling(d.x);
        }


        displaySequence(seq) {
            return this.width / seq > 5;
        }
        rectWidth(d) {
            return (this.scaling(d.y) - this.scaling(d.x));
        }

        rectX(object) {
            if (object.x === object.y) {
                return this.scaling(object.x-0.4);
            }
            return this.scaling(object.x);
        }

        rectWidth2(d){
            if (d.x === d.y) {
                if (this.scaling(d.x + 0.4) - this.scaling(d.x - 0.4) < 2) return 2;
                else return this.scaling(d.x + 0.4) - this.scaling(d.x - 0.4);
            }
            return (this.scaling(d.y) - this.scaling(d.x));
        }

        uniqueWidth(d) {
            return (this.scaling(1));
        }

        onFeatureSelected = function (listener) {
            this.svgElement.addEventListener(this.events.FEATURE_SELECTED_EVENT, listener);
        }
        onFeatureDeselected = function (listener) {
            this.svgElement.addEventListener(this.events.FEATURE_DESELECTED_EVENT, listener);
        }

        onZoom(listener) {
            this.svgElement.addEventListener(this.events.ZOOM_EVENT, listener);
        }

        addLevel(array) {
            var leveling = [];
            array.forEach(function (d) {
                if (leveling === []) {
                    leveling.push(d.y);
                    d.level = 0;
                } else {
                    var placed = false;
                    for (var k = 0; k < leveling.length; k++) {
                        if (d.x > leveling[k]) {
                            placed = true;
                            d.level = k;
                            leveling[k] = d.y;
                            break;
                        }
                    }
                    if (placed === false) {
                        leveling.push(d.y);
                        d.level = leveling.length - 1;
                    }
                }
            });
            return leveling.length;
        }

        addLevelToBond(array) {
            var leveling = [];
            var newArray = [];
            array.forEach(function (d) {
                if (leveling === []) {
                    leveling.push(d[2].x);
                    d[1].y = 1;
                } else {
                    var placed = false;
                    for (var k = 0; k < leveling.length; k++) {
                        if (d[0].x > leveling[k]) {
                            placed = true;
                            d[1].y = k + 1;
                            leveling[k] = d[2].x;
                            break;
                        }
                    }
                    if (placed === false) {
                        leveling.push(d[2].x);
                        d[1].y = leveling.length;
                    }
                }
            });
            return leveling.length;
        }

        sbcRip: any = null;
        
        shadeBlendConvert(p:number, from:string, to:string=null) {

            var self=this;

            if(typeof(p)!="number"||p<-1||p>1||typeof(from)!="string"||(from[0]!='r'&&from[0]!='#')||(typeof(to)!="string"&&typeof(to)!="undefined"))return null; //ErrorCheck
            if(!this.sbcRip)this.sbcRip=function(d){
                var l=d.length,RGB=new Object();
                if(l>9){
                    d=d.split(",");
                    if(d.length<3||d.length>4)return null;//ErrorCheck
                    RGB[0]=i(d[0].slice(4)),RGB[1]=i(d[1]),RGB[2]=i(d[2]),RGB[3]=d[3]?parseFloat(d[3]):-1;
                }else{
                    if(l==8||l==6||l<4)return null; //ErrorCheck
                    if(l<6)d="#"+d[1]+d[1]+d[2]+d[2]+d[3]+d[3]+(l>4?d[4]+""+d[4]:""); //3 digit
                    d=i(d.slice(1),16),RGB[0]=d>>16&255,RGB[1]=d>>8&255,RGB[2]=d&255,RGB[3]=l==9||l==5?r(((d>>24&255)/255)*10000)/10000:-1;
                }
                return RGB;}
            var i=parseInt,r=Math.round,h=from.length>9,h=typeof(to)=="string"?to.length>9?true:to=="c"?!h:false:h,b=p<0,p=b?p*-1:p,to=to&&to!="c"?to:b?"#000000":"#FFFFFF",f=self.sbcRip(from),t=self.sbcRip(to);
            if(!f||!t)return null; //ErrorCheck
            if(h)return "rgb("+r((t[0]-f[0])*p+f[0])+","+r((t[1]-f[1])*p+f[1])+","+r((t[2]-f[2])*p+f[2])+(f[3]<0&&t[3]<0?")":","+(f[3]>-1&&t[3]>-1?r(((t[3]-f[3])*p+f[3])*10000)/10000:t[3]<0?f[3]:t[3])+")");
            else return "#"+(0x100000000+(f[3]>-1&&t[3]>-1?r(((t[3]-f[3])*p+f[3])*255):t[3]>-1?r(t[3]*255):f[3]>-1?r(f[3]*255):255)*0x1000000+r((t[0]-f[0])*p+f[0])*0x10000+r((t[1]-f[1])*p+f[1])*0x100+r((t[2]-f[2])*p+f[2])).toString(16).slice(f[3]>-1||t[3]>-1?1:3);
        }

        addXAxis(position) {
            var self = this;

            this.svgContainer.append("g")
                .attr("class", "x axis Xaxis")
                .attr("transform", "translate(0," + (position + 20) + ")")
                .call(self.xAxis);
        }

        updateXaxis(position) {
            this.svgContainer.selectAll(".Xaxis")
                .attr("transform", "translate(0," + (position + 20) + ")")
        }

        updateSVGHeight(position) {
            this.svg.attr("height", position + 60 + "px");
            this.svg.select("clippath rect").attr("height", position + 60 + "px");
        }

        addYAxis() {
            this.yAxisSVG = this.svg.append("g")
                .attr("class", "pro axis")
                .attr("transform", "translate(0," + this.margin.top + ")");
            this.updateYaxis();
        }

        updateYaxis() {

            var self = this;

            this.yAxisSVGgroup = self.yAxisSVG
                .selectAll(".yaxis")
                .data(self.yData)
                .enter()
                .append("g");

            this.yAxisSVGgroup
                .append("polygon") // attach a polygon
                .attr("class", function (d) {
                    if (d.filter) return d.filter.split(" ").join("_") + "Arrow";
                    return "Arrow";
                })
                .style("stroke", "") // colour the line
                .style("fill", "#DFD5D3") // remove any fill colour
                .attr("points", function (d) {
                    return (self.margin.left - 105) + "," + (d.y - 3) + ", " + (self.margin.left - 105) + "," + (d.y + 12) + ", " + (self.margin.left - 15) + "," + (d.y + 12) + ", " + (self.margin.left - 7) + "," + (d.y + 4.5) + ", " + (self.margin.left - 15) + "," + (d.y -3); // x,y points
                });

            this.yAxisSVGgroup
                .append("text")
                .attr("class", "yaxis")
                .attr("text-anchor", "start")
                .attr("x", function () {
                    return self.margin.left - 102
                })
                .attr("y", function (d) {
                    return d.y + 8
                })
                .text(function (d) {
                    return d.title
                });
        }

        forcePropagation(item) {
            var self=this;


            item.on('mousedown', function () {
                var brush_elm = self.svg.select(".brush").node();

                var eventData = new MouseEventData();
                eventData.pageX = d3.event.pageX;
                eventData.clientX = d3.event.clientX;
                eventData.pageY = d3.event.pageY;
                eventData.clientY = d3.event.clientY;

                var new_click_event = new MouseEvent('mousedown', eventData);

                if (brush_elm) {
                    brush_elm.dispatchEvent(new_click_event);
                }
            });
        }

        /** export to new utils file  */

        showFilteredFeature(className, color, baseUrl){
            var featureSelected = this.yAxisSVG.selectAll("."+className+"Arrow");
            var minY = this.margin.left - 105;
            var maxY = this.margin.left - 7;

            var gradient = this.svg
                .append("linearGradient")
                .attr("y1", "0")
                .attr("y2", "0")
                .attr("x1", minY)
                .attr("x2", maxY)
                .attr("id", "gradient")
                .attr("spreadMethod", "pad")
                .attr("gradientUnits", "userSpaceOnUse");

            gradient
                .append("stop")
                .attr("offset", "0.3")
                .attr("stop-color", "#DFD5D3")
                .attr("stop-opacity", 1);


            var redgrad = gradient
                .append("stop")
                .attr("offset", "1")
                .attr("stop-opacity", 1)
                .attr("stop-color", "#DFD5D3");

            redgrad
                .attr("stop-color", color);

            var url_gradient = "url(#gradient)";
            var url_dropshadow = "url(#dropshadow)";
            if (baseUrl) {
                url_gradient = "url(" + baseUrl + "#gradient)";
                url_dropshadow = "url(" + baseUrl +"#dropshadow)";
            }

            var self=this;
            var selection = this.yAxisSVG.selectAll("."+className+"Arrow")
                .style("fill", url_gradient)
                .style("stroke", "")
                .attr("filter", url_dropshadow);
            selection
                .attr("points", function (d) {
                    return (self.margin.left - 105) + "," + (d.y - 3) + ", " + (self.margin.left - 105) + "," + (d.y + 12) + ", " + (self.margin.left - 10) + "," + (d.y + 12) + ", " + (self.margin.left - 2) + "," + (d.y + 4.5) + ", " + (self.margin.left - 10) + "," + (d.y -3); // x,y points
                });
        }


        hideFilteredFeature(className){
            var self=this;

            this.yAxisSVG.selectAll("."+className+"Arrow")
                .style("fill", "rgba(95,46,38,0.2)")
                .attr("filter", "")
                .attr("points", function (d) {
                    return (self.margin.left - 105) + "," + (d.y - 3) + ", " + (self.margin.left - 105) + "," + (d.y + 12) + ", " + (self.margin.left - 15) + "," + (d.y + 12) + ", " + (self.margin.left - 7) + "," + (d.y + 4.5) + ", " + (self.margin.left - 15) + "," + (d.y -3); // x,y points
                });
        };

        gBrush;

        addBrush() {
            var self=this;
            this.gBrush = this.svgContainer.append("g")
                .attr("class", "brush")
                .call(self.brush);

            this.gBrush.selectAll("rect")
                .attr('height', self.Yposition + 50);
        }
        
        zoom(start, end){

            var zoomInside = this.current_extend.start<start && this.current_extend.end>end;
            if (!zoomInside) {
                this.svgContainer.selectAll(".seqGroup").remove();
            }
            this.brush.extent([start,end]);
            this.brushend();
        }

        resetZoom(start, end){
            this.resetAll();
        }

        brushEmpty()
        {
            return d3.event.selection === null;
        }

        brushend() {

            var self = this;

            if (this.brushEmpty())
            {
                return;
            }


            if (!d3.event.sourceEvent) return; // Only transition after input.

            var eventType = d3.event.sourceEvent.type;

            console.log("brushend event function");
            console.log(eventType);

            //if (eventType != "brushend") return;

            console.log("After safety check");

            self.d3div.selectAll('div.selectedRect').remove();
            if (Object.keys(self.featureSelected).length !== 0 && self.featureSelected.constructor === Object) {
                d3.select(self.featureSelected.id).style("fill", self.featureSelected.originalColor);
                self.featureSelected = {};
                if (CustomEvent) {
                    var event = new CustomEvent(self.events.FEATURE_DESELECTED_EVENT, {
                        detail: {info:"feature-deselected"}
                    });
                    self.svgElement.dispatchEvent(event);
                } else {
                    console.warn("CustomEvent is not defined....");
                }
                if (self.trigger) self.trigger(self.events.FEATURE_DESELECTED_EVENT, {info:"feature-deselected"});
            }
            // Check if brush is big enough before zooming

            var extent = d3.event.selection;
            console.log(extent);
            var extentLength = 0;

            if (extent != null)
            {
                extent = [self.scalingPosition(extent[0]), self.scalingPosition(extent[1])]
                //var extent = self.brush.extent();
                extentLength = Math.abs(extent[0] - extent[1]);
            }

            console.log(extentLength);
            var brushEmpty = this.brushEmpty();
            console.log(brushEmpty);

            var seq = self.displaySequence(extentLength);
            if (!brushEmpty && extentLength > self.zoomMax) {

                console.log("Going into if");

                if (extent[0] < extent[1]) var start = extent[0] - 1,
                end = parseInt(extent[1] + 1);
                else var start = parseInt(extent[1] + 1),
                    end = extent[0] - 1;

                self.current_extend.length = extentLength;
                var zoomScale = (self.fvLength / extentLength).toFixed(1);
                $(self.viewerDiv).find(".zoomUnit").text(zoomScale.toString());
                
//                scaling.range([5,width-5]); 
                if (self.SVGOptions.showSequence && !(self.intLength) && seq && self.svgContainer.selectAll(".AA").empty()) {
                    self.current_extend = { 
                    length : extentLength,
                    start : start,
                    end : end
                    }
                    self.seqShift = start;
                    self.svgContainer.selectAll(".sequenceLine").remove();
                    self.fillSVG.sequence(self.sequence.substring(start-1, end), 20, self.seqShift-1);
                }

                //modify scale
//                scaling.range([5,width-5]);
                self.scaling.domain(extent);
                self.scalingPosition.range(extent);
                var currentShift = self.seqShift ? self.seqShift : self.offset.start;
                

                self.transition_data(self.features, currentShift);
                self.reset_axis();

                if (CustomEvent) {
                  self.svgElement.dispatchEvent(new CustomEvent(
                    self.events.ZOOM_EVENT,
                    {detail: { start: start, end: end, zoom: zoomScale }}
                    ));
                }
                if (self.trigger) self.trigger(self.events.ZOOM_EVENT, {
                            start: start,
                            end: end,
                            zoom: zoomScale
                        });

                //rectsPep2.classed("selected", false);
                //d3.select(self.div).selectAll(".brush").call(() => this.clearBrush());
            } else {
                //d3.select(self.div).selectAll(".brush").call(() => this.clearBrush());
                //resetAll();
            }

            this.clearBrush();
        }

        clearBrush(elem=null)
        {

            if (elem == null)
            {
                elem = this.gBrush;//d3.select('.brush').node() as HTMLElement;
            }
            console.log("clear brush");
            console.log(elem);
            this.gBrush.call(this.brush.move, null);



        }
//        
        trigger: any;

        resizeCallback(){
            this.updateWindow();
        }
                
        updateWindow(){
//            var new_width = $(div).width() - margin.left - margin.right - 17;
//            var width_larger = (width < new_width);
            var self=this;            

            self.width = $(self.viewerDiv).width() - self.margin.left - self.margin.right - 17;
            self.d3div.select("svg")
                .attr("width", self.width + self.margin.left + self.margin.right);
            self.d3div.select("#clip>rect").attr("width", self.width);
            if (self.SVGOptions.brushActive) {
                self.d3div.select(".selection").attr("width", self.width);
            }
            self.d3div.selectAll(".brush").call(() => this.clearBrush());
            
            //var currentSeqLength = svgContainer.selectAll(".AA").size();
            var seq = self.displaySequence(self.current_extend.length);
            if (self.SVGOptions.showSequence && !(self.intLength)){
                if (seq === false && !self.svgContainer.selectAll(".AA").empty()) {
                    self.svgContainer.selectAll(".seqGroup").remove();
                    self.fillSVG.sequenceLine();
                }
                else if (seq === true && self.svgContainer.selectAll(".AA").empty()){
                    self.svgContainer.selectAll(".sequenceLine").remove();
                    self.fillSVG.sequence(self.sequence.substring(self.current_extend.start-1, self.current_extend.end), 20, self.current_extend.start-1);
                    
                }
            }
            
            self.scaling.range([5,self.width-5]);
            self.scalingPosition.domain([0, self.width]);
            
            self.transition_data(self.features, self.current_extend.start);
            self.reset_axis();
            
        }

        // If brush is too small, reset view as origin
        resetAll() {

            //reset scale
            var self=this;

            $(".zoomUnit").text("1");
            self.scaling.domain([self.offset.start, self.offset.end]);
            self.scalingPosition.range([self.offset.start, self.offset.end]);
            var seq = self.displaySequence(self.offset.end - self.offset.start);
            
            if (self.SVGOptions.showSequence && !(self.intLength)){
                if (seq === false && !self.svgContainer.selectAll(".AA").empty()){
                    self.svgContainer.selectAll(".seqGroup").remove();
                    self.fillSVG.sequenceLine();
                }
                else if (self.current_extend.length !== self.fvLength && seq === true && !self.svgContainer.selectAll(".AA").empty()) {
                    self.svgContainer.selectAll(".seqGroup").remove();
                    self.fillSVG.sequence(self.sequence.substring(self.offset.start-1,self.offset.end), 20, self.offset.start);
                }
            }

            self.current_extend={ 
                    length : self.offset.end-self.offset.start,
                    start : self.offset.start,
                    end : self.offset.end
                };
            self.seqShift=0;
            
            self.transition_data(self.features, self.offset.start);
            self.reset_axis();

            // Fire Event
            if (CustomEvent) {
                self.svgElement.dispatchEvent(new CustomEvent(self.events.ZOOM_EVENT,
                { detail: { start: 1, end: self.sequence.length, zoom: 1 }}));
            };
            if (self.trigger) self.trigger(self.events.ZOOM_EVENT, {
                            start: 1,
                            end: self.sequence.length,
                            zoom: 1
                        });

            self.d3div.selectAll(".brush").call(() => this.clearBrush());
        }

        /** export to new axis file? */
        reset_axis() {
            var self=this;
            self.svgContainer
                .transition().duration(500)
                .select(".x.axis")
                .call(self.xAxis);
        }

        addVerticalLine() {
            var self=this;
            var vertical = d3.select(".chart")
                .append("div")
                .attr("class", "Vline")
                .style("position", "absolute")
                .style("z-index", "19")
                .style("width", "1px")
                .style("height", (self.Yposition + 50) + "px")
                .style("top", "30px")
                // .style("left", "0px")
                .style("background", "#000");

            var mouseHandler = function () {
                var mousex = d3.mouse(this)[0] - 2;
                vertical.style("left", mousex + "px")
            };

            d3.select(".chart")
                .on("mousemove.Vline", mouseHandler);
            //.on("click", function(){
            //    mousex = d3.mouse(this);
            //    mousex = mousex[0] + 5;
            //    vertical.style("left", mousex + "px")});
        }

        addRectSelection(svgId) {
            var self=this;
            var div = self.d3div;

            var featSelection = d3.select(svgId);
            var elemSelected: Array<any> = featSelection.data();
            var xTemp;
            var yTemp;
            var xRect;
            var widthRect;
            var svgWidth = self.SVGOptions.brushActive ? d3.select(".selection").attr("width") : self.svgContainer.node().getBBox().width;

            // ATTENTION TODO THIS MIGHT INTERFERE WITH OTHER FV!
            //d3.select('body').selectAll('div.selectedRect').remove();
            self.d3div.selectAll('div.selectedRect').remove();

        
            var objectSelected = {type:featSelection[0][0].tagName, color:featSelection.style("fill")};
            self.tooltip.colorSelectedFeat(svgId, objectSelected);

            // Append tooltip
            var selectedRect = div.append('div')
                .attr('class', 'selectedRect');

            if (elemSelected[0].length === 3) {
                xTemp = elemSelected[0][0].x;
                yTemp = elemSelected[0][1].x;
            } else if (elemSelected[0].x === elemSelected[0].y) {
                xTemp = elemSelected[0].x - 0.5;
                yTemp = elemSelected[0].y + 0.5;
            } else {
                xTemp = elemSelected[0].x;
                yTemp = elemSelected[0].y;
            }
            if (self.scaling(xTemp) < 0) {
                xRect = self.margin.left;
                widthRect = (self.scaling(yTemp));
            } else if (self.scaling(yTemp) > svgWidth) {
                xRect = self.scaling(xTemp) + self.margin.left;
                widthRect = svgWidth - self.scaling(xTemp);
            } else {
                xRect = self.scaling(xTemp) +self. margin.left;
                widthRect = (self.scaling(yTemp) - self.scaling(xTemp));
            }

            selectedRect.style("left", xRect + "px");
            selectedRect.style("top", 60 + "px");
            selectedRect.style("width", widthRect + "px");
            selectedRect.style("height", (self.Yposition + 50) + "px");
            selectedRect.style("position", "absolute");
            selectedRect.style("z-index", xRect + "px");
            selectedRect.style("box-shadow", "0 1px 2px 0 #656565");
            selectedRect.style("background-color", 'rgba(0, 0, 0, 0.2)');
        }

        initSVG(div, options: 
            {showAxis: boolean,
                showSequence:boolean,
                brushActive:boolean,
                verticalLine:boolean,
                toolbar: boolean,
                bubbleHelp:boolean,
                unit: string,
                zoomMax: number,
                animation:boolean,
                positionWithoutLetter: boolean,
                dottedSequence: boolean
            }) {

            var self=this;

            if (typeof options === 'undefined') {
                var options = {
                    'showAxis': false,
                    'showSequence': false,
                    'brushActive': false,
                    'verticalLine': false,
                    'toolbar': false,
                    'bubbleHelp': false,
                    'unit': "units",
                    'zoomMax': 50,
                    'animation': false,
                    'dottedSequence': true,
                    'positionWithoutLetter': false
                }
            }

            if (!$.fn.popover) {
                options.bubbleHelp = false;
                console.warn("The bubble help requires tooltip and popover bootstrap js libraries. The feature viewer will continue to work, but without the info bubble");
            }

            // Create SVG
            if (options.zoomMax) {
                self.zoomMax = options.zoomMax;
            }
            if (!options.unit) {
                options.unit = "units";
            }
            if (options.animation) {
                self.animation = options.animation;
            }

            if (options.toolbar === true) {
                
                var headerOptions = $(div).find(".svgHeader").length ? d3.select(div).select(".svgHeader") : d3.select(div).append("div").attr("class", "svgHeader");
                
//                if (options.toolbarTemplate && options.toolbarTemplate === 2) {

                    if (!$(div).find('.header-position').length) {
                        var headerPosition = headerOptions
                            .append("div")
                            .attr("class", "header-position")
                            .style("display", "inline-block")
                            .style("margin", "15px 10px 0px")
                            .style("padding", "0px")
                            .style("line-height","32px");
                        headerPosition
                            .append("div")
                            .attr("class", "position-label")
                            .style("padding", "0px 5px")
                            .style("display", "inline-block")
                            .style("padding", "0px")
                            .style("font-weight","700")
                            .text("Position  :  ");
                        headerPosition
                            .append("div")
                            .style("display", "inline-block")
                            .style("padding", "0px")
                            .style("padding-left", "5px")
                            .append("div")
                            .style("min-width","50px")
                            .attr("id", "zoomPosition")
                            .text("0");
                    }
                    if (!$(div).find(' .header-zoom').length) {
                        var headerZoom = headerOptions
                            .append("div")
                            .attr("class", "header-zoom")
                            .style("display", "inline-block")
                            .style("margin", "15px 0px 0px")
                            .style("padding", "0px")
                            .style("line-height","32px");
                        headerZoom
                            .append("div")
                            .attr("class", "zoom-label")
                            .style("padding", "0px 5px")
                            .style("display", "inline-block")
                            .style("padding", "0px")
                            .style("font-weight","700")
                            .text("Zoom : ");

                        headerZoom
                            .append("div")
                            .style("display", "inline-block")
                            .style("padding", "0px")
                            .append("div")
                            .style("min-width","50px")
                            .style("padding-left", "5px")
                            .append("span")
                            .text("x ")
                            .append("span")
                            .attr("class", "zoomUnit")
                            .text("1");
                    }
//                }
//                else{
//                    if (!$(div + ' .header-zoom').length) {
//                        var headerZoom = headerOptions
//                            .append("div")
//                            .attr("class", "panel panel-default header-zoom")
//                            .style("display", "inline-block")
//                            .style("width", "150px")
//                            .style("margin", "20px 0px 0px")
//                            .style("padding", "0px");
//                        headerZoom
//                            .append("div")
//                            .attr("class", "panel-heading")
//                            .style("padding", "0px 15px")
//                            .style("border-right", "1px solid #DDD")
//                            .style("display", "inline-block")
//                            .style("width", "80px")
//                            .append("h5")
//                            .style("padding", "0px")
//                            .style("height", "10px")
//                            .style("color", "#777")
//                            .text("ZOOM");
//                        headerZoom
//                            .append("div")
//                            .attr("class", "panel-body")
//                            .style("display", "inline-block")
//                            .style("padding", "0px")
//                            .append("h5")
//                            .style("padding-left", "15px")
//                            .style("height", "10px")
//                            .text("x ")
//                            .append("span")
//                            .attr("class", "zoomUnit")
//                            .text("1");
//                    }
//                    if (!$(div + ' .header-position').length) {
//                        var headerPosition = headerOptions
//                            .append("div")
//                            .attr("class", "panel panel-default header-position")
//                            .style("display", "inline-block")
//                            .style("width", "175px")
//                            .style("margin", "20px 20px 0px")
//                            .style("padding", "0px");
//                        headerPosition
//                            .append("div")
//                            .attr("class", "panel-heading")
//                            .style("padding", "0px 15px")
//                            .style("border-right", "1px solid #DDD")
//                            .style("display", "inline-block")
//                            .append("h5")
//                            .style("padding", "0px")
//                            .style("height", "10px")
//                            .style("color", "#777")
//                            .text("POSITION");
//                        headerPosition
//                            .append("div")
//                            .attr("class", "panel-body")
//                            .style("display", "inline-block")
//                            .style("padding", "0px")
//                            .append("h5")
//                            .style("padding-left", "15px")
//                            .style("height", "10px")
//                            .append("span")
//                            .attr("id", "zoomPosition")
//                            .text("0");
//                    }
//                }

                var headerZoom = d3.select(div).select(' .header-zoom');
                console.log("Header Zoom Element");
                console.log(headerZoom);
                console.log(headerOptions);

                var headerZoom = $(div).find('.header-zoom').length ? d3.select(div).select('.header-zoom') : headerOptions;
                if (options.bubbleHelp === true) {
                    if (!$(div).find('.header-help').length) {
                        var helpContent = "<div><strong>To zoom in :</strong> Left click to select area of interest</div>" +
                            "<div><strong>To zoom out :</strong> Right click to reset the scale</div>" +
                            "<div><strong>Zoom max  :</strong> Limited to <strong>" + self.zoomMax.toString() + " " + options.unit +"</strong></div>";
//                        var headerHelp = headerOptions
                        var headerHelp = headerZoom
                            .append("div")
//                            .insert("div",":first-child")
//                            .attr("class", "pull-right")
                            .style("display", "inline-block")
//                            .style("margin", "15px 35px 0px 0px")
                            .style("margin", "0px")
                            .style("margin-right", "5px")
//                            .style("line-height","32px")
                            .style("padding", "0px");
                        var buttonHelp = headerHelp
                            .append("a")
                            .attr("type", "button")
                            .attr("class", "header-help")
                            .attr("data-toggle", "popover")
                            .attr("data-placement", "auto left")
                            .attr("title", "Help")
                            //.attr("data-content", helpContent)
                            .style("font-size", "14px");
//                            .style("margin-bottom", "2px");
                        buttonHelp
                            .append("span")
                            .attr("class", "label label-as-badge label-info")
                            .style("font-weight","500")
//                            .style("border-radius","3px")
                            .style("border-radius","3px")
//                            .style("background-color","#f8f8f8")
//                            .style("background-color","#108D9F")
//                            .style("border","1px solid #ddd")
//                            .style("border","1px solid #0C6B78")
//                            .style("color","#777")
                            .style("color", "#000000")
                            .style("box-shadow","inset 0px 0px 4px rgba(0,0,0,0.10)")
                            //.style("color","#fff")
//                            .style("padding","2px 6px")
                            .html("<span class='state'>Show</span> help");
                        $(function () {
                            //$('[data-toggle="popover"]').popover({html: true});
                            //var baseContID = self.div.substring(1);
                            //console.log("base div id: " + baseContID);

                            var helpButton = $(div).find('.header-help');

                            helpButton.popover({
                                container: "body",
                                html: true,
                                placement: 'right',
                                content: function () {
                                  return helpContent;
                                }
                              });

                              helpButton.on('hide.bs.popover', function () {
                                helpButton.find(".state").text("Show");
                            });
                            helpButton.on('show.bs.popover', function () {
                                helpButton.find(".state").text("Hide");
                            });
                        })
                    }
                }
            }
            
            self.svg = d3.select(div).append("svg")
                .attr("width", self.width + self.margin.left + self.margin.right)
                .attr("height", self.height + self.margin.top + self.margin.bottom)
                .style("z-index", "2")
                .on("contextmenu", function (d, i) {
                    d3.event.preventDefault();
                    self.resetAll();
                    // react on right-clicking
                });
            self.svgElement = self.viewerDiv.getElementsByTagName("svg")[0];

            self.svgContainer = self.svg
                .append("g")
                .attr("transform", "translate(" + self.margin.left + "," + self.margin.top + ")");

            //Create Clip-Path
            var defs = self.svgContainer.append("defs");

            defs.append("clipPath")
                .attr("id", "clip")
                .append("rect")
                .attr("width", self.width)
                .attr("height", self.height);

            var filter = defs.append("filter")
                .attr("id", "dropshadow")
                .attr("height", "200%");

            filter.append("feGaussianBlur")
                .attr("in", "SourceAlpha")
                .attr("stdDeviation", 3)
                .attr("result", "blur");
            filter.append("feOffset")
                .attr("in", "blur")
                .attr("dx", -2)
                .attr("dy", 2)
                .attr("result", "offsetBlur");

            var feMerge = filter.append("feMerge");

            feMerge.append("feMergeNode")
                .attr("in", "offsetBlur");
            feMerge.append("feMergeNode")
                .attr("in", "SourceGraphic");

                self.svgContainer.on('mousemove', function () {

                var bgnodes = d3.selectAll(".selection").nodes();
                var elem: HTMLElement = bgnodes[0] as HTMLElement;

                //console.log("initSVG");
                //console.log(bgnodes);
                //console.log(elem);
                //console.log(self.SVGOptions);

                var absoluteMousePos = self.SVGOptions.brushActive ? d3.mouse(elem) : d3.mouse(self.svgContainer);       
                var pos = Math.round(self.scalingPosition(absoluteMousePos[0]));
                if (!options.positionWithoutLetter) {
                    pos += self.sequence[pos-1] || "";
                }
                $(div).find("#zoomPosition").text(pos);
            });
            
            if (typeof options.dottedSequence !== "undefined"){
                self.SVGOptions.dottedSequence = options.dottedSequence;
            }
            if (options.showSequence && !(self.intLength)) {
                self.SVGOptions.showSequence = true;
                if (self.displaySequence(self.offset.end - self.offset.start)) {
                    self.fillSVG.sequence(self.sequence.substring(self.offset.start-1, self.offset.end), self.Yposition, self.offset.start);
                }
                else{
                    self.fillSVG.sequenceLine();
                }
                self.features.push({
                    data: self.sequence,
                    name: "Sequence",
                    className: "AA",
                    color: "black",
                    type: "text"
                });
                self.yData.push({
                    title: "Sequence",
                    y: self.Yposition - 8
                });
            }
            if (options.showAxis) 
            {
                self.addXAxis(self.Yposition);
            }

            self.addYAxis();

            if (options.brushActive) {
                self.SVGOptions.brushActive = true;
                self.enableZoom = true;
                self.addBrush();
            }
            if (options.verticalLine) {
                self.SVGOptions.verticalLine = true;
                self.addVerticalLine();
            }

            self.updateSVGHeight(self.Yposition);

        }

        addFeature(object) {
            this.Yposition += 20;
            this.features.push(object);
            this.fillSVG.typeIdentifier(object);
            this.updateYaxis();
            this.updateXaxis(this.Yposition);
            this.updateSVGHeight(this.Yposition);
            if (this.SVGOptions.brushActive) {
                this.svgContainer.selectAll(".brush rect")
                    .attr('height', this.Yposition + 50);
            }
            if (this.SVGOptions.verticalLine) d3.selectAll(".Vline").style("height", (this.Yposition + 50) + "px");

            var allElems = d3.selectAll(".element").nodes();
            console.log("Add Feature");
            console.log(object);
            console.log(allElems);

            if ((allElems) && (allElems.length > 1500)) this.animation = false;

        }
        
        clearInstance (){
            $(window).off("resize", this.resizeCallback.bind(this));
            this.svg = null;
            this.svgElement = null;
            this.svgContainer = null;
            this.yAxisSVGgroup = null;
            this.yAxisSVG = null;
            this.features = null;
            this.sbcRip = null;

        }

}