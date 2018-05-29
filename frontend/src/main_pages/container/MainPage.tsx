import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';

import D3ParallelLinesGraph from '../components/D3ForceParallelLineViewer';


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
     */
    render() {

        var graphData = {
            nodes: [
                {id: 'CXCR4', label: 'CXCR4'},
                {id: 'CCL5', label: 'CCL5'}
            ],
            links: [
                {source: 0, target: 1, evidence: 20, predicted: 20}
            ]
        }


        return (
            <div>
                <span>
                    <h4>General Information</h4>
                    <ul>
                        <li>Number of Genes</li>
                        <li>Number of miRNA</li>
                        <li>Number of lncRNAs</li>
                        <li>Number if Gene/miRNA Interactions</li>
                        <li>Number of lncRNA/miRNA Interactions</li>
                    </ul>
                </span>
                <span>
                    <h4>Included Databases</h4>
                    <ul>
                        <li>mirTarBase</li>
                        <li>miRecords</li>
                        <li>PubMeds</li>
                    </ul>
                </span>
                <D3ParallelLinesGraph graph={graphData} id="d3front"/>
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
                    title="Welcome!"
                    subtitle="miRExplore v0.1"
                />
                <CardText >
                    <p>miRExplore v.01 online.</p>
                    <p>Explore the database in the explore tab</p>
                    <p>Provide some statistics of the project here</p>

                </CardText>
            </Card>
                        <Card>
                        <CardTitle
                            title="miRExplore"
                            subtitle="Statistics"
                        />
                        <CardText >
                            <MainStatus/>
                        </CardText>
                    </Card>
            </div>

        );
    }

}