import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';

import D3SVGParallelLinesGraph from '../components/D3SVGForceParallelLines';
import MEFeatureViewer from '../components/MEFeatureViewer';
import SankeyChart from '../components/SankeyChart';

export class MainStatus extends React.Component<{},{}> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

    }

    /**
     * This method will be executed after initial rendering.
     */
    componentDidMount() {

    }

    /**
     * Render the component.
     * <MEFeatureViewer features={[]} location={[]}/>
     * 
     *                 <MEFeatureViewer features={[]} location={[]}/>
                <D3SVGParallelLinesGraph graph={graphData} id="d3front"/>
     */
    render() {

        var graphData = {
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


        return (
            <div>
                <span>
                    <h4>General Information</h4>
                    <p>
                        neutrophil TM 
                    </p>
                </span>
                <span>
                    <h4>Included Databases</h4>
                    <ul>
                        <li><b>Pubmed</b></li>
                    </ul>
                    <h4>To be done</h4>
                    <ul>
                        <li><b>PMC</b>: fulltexts and TM by section</li>
                    </ul>
                </span>

                <SankeyChart graph={graphData} id="d3front"/>
                

            </div>
        );
    }

}


export class WelcomePage extends React.Component<{},{}> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

    }

    /**
     * This method will be executed after initial rendering.
     */
    componentDidMount() {

    }

    /**
     * Render the component.
     */
    render() {


        return (

            <div>
            <Card style={{marginBottom: "20px"}}>
                <CardTitle
                    title="neutrophilTM"
                    subtitle="v1.0"
                />
                <CardText >
                    <p>neutrophil TM</p>
                    <p>Welcome!</p>
                    <p>Explore the database in the explore tab</p>
                </CardText>
            </Card>
                    <Card>
                        <CardText >
                            <MainStatus/>
                        </CardText>
                    </Card>
            </div>

        );
    }

}