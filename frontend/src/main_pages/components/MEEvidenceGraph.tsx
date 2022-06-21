import * as React from "react"; 
import * as ReactDOM from "react-dom";
import FeatureViewer from './FeatureViewer';

import * as pyc from "pycollections";
import D3SVGParallelLinesGraph from '../components/D3SVGForceParallelLines';
import axios from 'axios';
import config from '../config';

export interface MEEvidenceGraphProps { data: any}
export interface MEEvidenceGraphState { graphData:any, graphInfo:any, graphCoExpression:any}



export default class MEEvidenceGraph extends React.Component<MEEvidenceGraphProps, MEEvidenceGraphState> {

    constructor(props) {
        super(props);

        console.log(this.props);



    }

    readonly state = {
        'graphData': {
            'nodes': [],
            'links': []
        },
        'graphInfo':{},
        'graphCoExpression': []
    };

    componentDidUpdate(prevProps, prevState, snapshot) {

        var self=this;

        if (prevProps != this.props)
        {

            this.calcNewGraphData();
            this.getNewExpression();

        }


    }

    componentDidMount()
    {
        this.calcNewGraphData();
        this.getNewExpression();

    }

    calcNewGraphData()
    {
        console.log("ME Evidence Graph Render");

        var graphData = {
            nodes: [],
            links: []
        };

        console.log(this.props.data);
        console.log(this.props.data.length);


        for (let i = 0; i < this.props.data.length; ++i)
        {
            var dataElem = this.props.data[i];
            var source = dataElem.lid;
            var target = dataElem.rid;

            var sourceIdx = -1;
            var targetIdx = -1;

            /**
             * 
             * FIND SOURCE
             */
            var foundSrcElems = graphData.nodes.filter((elem) => elem.id == source)

            if (foundSrcElems.length == 0)
            {
                // add node
                let nodeElem = {id: dataElem.lid, group: dataElem.ltype, name: dataElem.lid, idx: graphData.nodes.length};
                graphData.nodes.push(nodeElem);
                foundSrcElems.push(nodeElem);
            }

            sourceIdx = graphData.nodes.indexOf(foundSrcElems[0]);

            /**
             * 
             * FIND TARGET
             */

            var foundTgtElems = graphData.nodes.filter((elem) => elem.id == target)

            if (foundTgtElems.length == 0)
            {
                // add node
                let nodeElem = {id: dataElem.rid, group: dataElem.rtype, name: dataElem.rid, idx: graphData.nodes.length};
                graphData.nodes.push(nodeElem);
                foundTgtElems.push(nodeElem);
            }

            targetIdx = graphData.nodes.indexOf(foundTgtElems[0]);

            if ((sourceIdx == -1) || (targetIdx == -1))
            {
                console.log("srcId or tgtIdx == -1");
                continue;
            }

            /**
             * 
             * PREPARE EDGE
             */

            var evidenceEdges = 0;
            var predictedEdges = 0;

            for (var j =0; j < dataElem.evidences.length; ++j)
            {
                var ev = dataElem.evidences[j];

                if (ev['data_source'] in ['pmid', ])
                {
                    evidenceEdges += 1;
                } else {
                    predictedEdges += 1;
                }
            }

            graphData.links.push({source: sourceIdx, target: targetIdx, group1: evidenceEdges, group2: predictedEdges, group3: 0})

        }

        if (this.state.graphCoExpression)
        {

            for (var i = 0; i < this.state.graphCoExpression.length; ++i)
            {
                var dataElem = this.state.graphCoExpression[i];
                var source = dataElem.source;
                var target = dataElem.target;
    
                var sourceIdx = -1;
                var targetIdx = -1;
    
                console.log(dataElem);
    
                /**
                 * 
                 * FIND SOURCE
                 */
                var foundSrcElems = graphData.nodes.filter((elem) => elem.id == source)
    
                if (foundSrcElems.length == 0)
                {
                    // add node
                    let nodeElem = {id: dataElem.source, name: dataElem.source, idx: graphData.nodes.length};
                    graphData.nodes.push(nodeElem);
                    foundSrcElems.push(nodeElem);
                }
    
                sourceIdx = graphData.nodes.indexOf(foundSrcElems[0]);
    
                /**
                 * 
                 * FIND TARGET
                 */
    
                var foundTgtElems = graphData.nodes.filter((elem) => elem.id == target)
    
                if (foundTgtElems.length == 0)
                {
                    // add node
                    let nodeElem = {id: dataElem.target, name: dataElem.target, idx: graphData.nodes.length};
                    graphData.nodes.push(nodeElem);
                    foundTgtElems.push(nodeElem);
                }
    
                targetIdx = graphData.nodes.indexOf(foundTgtElems[0]);
    
                if ((sourceIdx == -1) || (targetIdx == -1))
                {
                    console.log("srcId or tgtIdx == -1");
                    continue;
                }
    
                /**
                 * 
                 * PREPARE EDGE
                 */
    
                var evidenceEdges = 0;
                var predictedEdges = 0;

                var ecolor = "#CA8D43";

                var hasCoexpressionScore = false;

                Object.keys(dataElem).forEach( deKey => {
                    var result = deKey.match(/coexpression/i);

                    if (result)
                    {
                        hasCoexpressionScore = true;
                    }
                });

                if (hasCoexpressionScore)
                {
                    ecolor = "#FF6D43";
                }
        
                graphData.links.push({source: sourceIdx, target: targetIdx, group1: 0, group2: 0, group3: 1, group3color: ecolor})
            }

        }

        var allGeneNames = [];

        graphData.nodes.forEach(element => {
            allGeneNames.push(element.id);
        })

        var data = {
            "genes": allGeneNames
        }
        var self=this;

        axios.post(config.getRestAddress() + "/get_expression", data, config.axiosConfig)
        .then(function (response) {

            console.log("ME Interaction Graph: expression")
            console.log(response.data)

            var newGraphInfo = response.data;
            self.setState({graphInfo: newGraphInfo});                    
        })
        .catch(function (error) {
        console.log(error);
        self.setState({graphCoExpression: {}});
        });


        console.log("Prepared graphData")

        this.setState({"graphData": graphData});
    }

    getNewExpression()
    {

        console.log("Attempting to fetch expression");
        var data = {
            "genes": []
        };

        for (let i = 0; i < this.props.data.length; ++i)
        {
            var dataElem = this.props.data[i];
            var source = dataElem.lid;
            var target = dataElem.rid;


            if (data.genes.indexOf(source) < 0)
            {
                data.genes.push(source);
            }

            if (data.genes.indexOf(target) < 0)
            {
                data.genes.push(target);
            }

        }


        console.log(data);

        var self=this;
        axios.post(config.getRestAddress() + "/get_expression", data, config.axiosConfig)
        .then(function (response) {

            console.log("ME Interaction Graph: expression")
            console.log(response.data)

            var newGraphInfo = response.data;

            self.setState({graphInfo: newGraphInfo});

            var additionalGenes = [];

            Object.keys(newGraphInfo).forEach(element => {
                
                var geneInfos = newGraphInfo[element];

                geneInfos.forEach(geneInfo => {
                    
                    if ("left_gene" in geneInfo)
                    {
                        additionalGenes.push(geneInfo["left_gene"]);
                    }

                    if ("right_gene" in geneInfo)
                    {
                        additionalGenes.push(geneInfo["right_gene"]);
                    }

                    if ("overlap" in geneInfo)
                    {
                        additionalGenes = additionalGenes.concat(geneInfo["overlap"]);
                    }

                });


            });



            data.genes = data.genes.concat(additionalGenes);

            console.log("Coexpression Genes")            
            console.log(data)
            // here we need to add the LNC from the graphInfo
            axios.post(config.getRestAddress() + "/get_coexpression", data, config.axiosConfig)
            .then(function (response) {

                console.log("ME Interaction Graph: coexpression")
                console.log(response.data)

                self.setState({graphCoExpression: response.data});
                self.calcNewGraphData();
                
            })
            .catch(function (error) {
            console.log(error);
            self.setState({graphCoExpression: {}});
            });

            
            
        })
        .catch(function (error) {
          console.log(error);
          self.setState({graphInfo: {}});
        });

    }



    render() {       

        return (<div>
            <D3SVGParallelLinesGraph graph={this.state.graphData} graphInfo={this.state.graphInfo}/>
        </div>);
    }

}