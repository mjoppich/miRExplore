import * as React from "react"; 
import * as ReactDOM from "react-dom";
import FeatureViewer from './FeatureViewer';

import * as pyc from "pycollections";
import D3SVGParallelLinesGraph from '../components/D3SVGForceParallelLines';


export interface MEEvidenceGraphProps { data: any}
export interface MEEvidenceGraphState { newData: any, graphData:any}



export default class MEEvidenceGraph extends React.Component<MEEvidenceGraphProps, MEEvidenceGraphState> {

    constructor(props) {
        super(props);

        console.log(this.props);



    }

    componentWillMount()
    {
        this.setState({newData: [], graphData: []});
    }

    componentDidMount()
    {
        if (this.props.data)
        {
            this.setState({newData: this.props.data});
        }
    }

    componentWillReceiveProps(newprops)
    {
        console.log("ME EG will receive props");
        console.log(newprops);

        if (newprops.data)
        {
            this.setState({newData: newprops.data});
        }
    }



    render() {
        /**
         * 
         *         var graphData = {
            nodes: [
                {id: 'CXCR4', label: 'CXCR4'},
                {id: 'CX3CL1', label: 'CX3CL1'},
                {id: 'CCL5', label: 'CCL5'}
            ],
            links: [
                {source: 0, target: 1, evidence: 20, predicted: 50},
                {source: 1, target: 2, evidence: 40, predicted: 80},
                {source: 0, target: 2, evidence: 80, predicted: 90}
            ]
        }
         * 
         */

        console.log("ME Evidence Graph Render");

        var graphData = {
            nodes: [],
            links: []
        };

        console.log(this.state.newData);
        console.log(this.state.newData.length);


        for (let i = 0; i < this.state.newData.length; ++i)
        {
            var dataElem = this.state.newData[i];
            var source = dataElem.lid;
            var target = dataElem.rid;

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
                let nodeElem = {id: dataElem.lid, group: dataElem.ltype, label: dataElem.lid};
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
                let nodeElem = {id: dataElem.rid, group: dataElem.rtype, label: dataElem.rid};
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

            graphData.links.push({source: sourceIdx, target: targetIdx, evidence: evidenceEdges, predicted: predictedEdges})

        }

        console.log("Prepared graphData")
        console.log(graphData);

        return (<div>
            <D3SVGParallelLinesGraph graph={graphData} id="d3front"/>
        </div>);
    }

}