import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';

import D3SVGParallelLinesGraph from '../components/D3SVGForceParallelLines';
import MEFeatureViewer from '../components/MEFeatureViewer';

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
                        mingleRNA is an interactive web tool, which allows the user to browse miRNA - lncRNA - mRNA interactions in human and mouse. 
                    </p>
                </span>
                <span>
                    <h4>Included Databases</h4>
                    <ul>
                        <li><b>miRBase</b>: 1918 (human), 1227 (mouse)</li>
                        <li><b>GENCODE</b>: 1882 (human miRNA), 28469 (human lncRNA), 2203 (mouse miRNA), 17515 (mouse lncRNA)</li>
                        <li><b>NONCODE</b>: 172218 (human), 131699 (mouse)</li>
                        <li><b>Lncipedia</b>: 119439 (human), no mouse data</li>
                        <li><b>lncRNAdb</b>: 164 (human), 98 (mouse)</li>
                        <li><b>miRanda</b>:
                            <ul> 
                                <li>1097064 (miRNA - gene, human, all), 819078 (miRNA - gene, mouse, all)</li>
                                <li>~4 mil. (miRNA - lncRNA, human, filtered), ~1 mil. (miRNA - lncRNA, mouse, filtered)</li>
                            </ul>
                        </li>
                        <li><b>mirTarBase</b></li>
                        <li><b>miRecords</b></li>
                        <li><b>PubMeds</b></li>
                    </ul>
                </span>
                
                <MEFeatureViewer features={[]} location={[]}/>
                <D3SVGParallelLinesGraph graph={graphData} id="d3front"/>
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
                    title="mingleRNA"
                    subtitle="v0.1"
                />
                <CardText >
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