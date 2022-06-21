import * as React from "react"; 
import * as ReactDOM from "react-dom";
import FeatureViewer from './FeatureViewer';
import axios from 'axios';
import config from '../config';

import D3SVGParallelLinesGraph from './D3SVGForceParallelLines';

export interface NIInteractionNetworkProps { obolevels: any, data: any, graphtitle:string}
export interface NIInteractionNetworkState { newData: any, graphData:any}



export default class NISankeyChart extends React.Component<NIInteractionNetworkProps, NIInteractionNetworkState> {

    emptyGraph = {
        'nodes': [],
        'links': []
    };

    readonly state = {
        'graphData': this.emptyGraph,
        'newData': []
    };



    constructor(props) {
        super(props);

        console.log(this.props);

    }

    componentWillMount()
    {
    }

    componentDidMount()
    {
        console.log("NI iNetwork did mount");
        this.reloadGraphData();
    }


    componentDidUpdate(prevProps, prevState, snapshot) {

        if (prevProps != this.props)
        {
            this.reloadGraphData();
        }
    }

    reloadGraphData()
    {

        var sendData = this.props.data;

        if (sendData == null)
        {
            return;
        }

        var self=this;

        console.log("NI iNetwork reload graph data");
        console.log(sendData);

        sendData['obolevel'] = this.props.obolevels['cells'];
        sendData['messenger_obolevel'] = this.props.obolevels['messenger'];


        axios.post(config.getRestAddress() + "/interaction_network",sendData, config.axiosConfig)
        .then(function (response) {

            console.log("NI iNetwork data fetch")
            console.log(response.data)

            self.setState({graphData: response.data});
            
        })
        .catch(function (error) {
          console.log(error);
          self.setState({graphData: self.emptyGraph});
        });
    }

    render() {

        return (<div>
            <D3SVGParallelLinesGraph graph={this.state.graphData} graphtitle={this.props.graphtitle}/>
        </div>);
    }

}