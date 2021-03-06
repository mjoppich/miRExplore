import * as React from "react"; 
import * as ReactDOM from "react-dom";
import FeatureViewer from './FeatureViewer';
import axios from 'axios';
import config from '../config';
import D3SankeyChart from "./SankeyChart";


export interface NISankeyChartProps { obolevels: any, data: any, graphtitle:string}
export interface NISankeyChartState { newData: any, graphData:any}



export default class NISankeyChart extends React.Component<NISankeyChartProps, NISankeyChartState> {

    emptyGraph = {
        'nodes': [],
        'links': []
    };


    constructor(props) {
        super(props);

        console.log(this.props);

    }


    readonly state = {
        'graphData': this.emptyGraph,
        'newData': []
    };

    componentDidMount()
    {
        console.log("NI Sankey did mount");
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

        sendData['obolevel'] = this.props.obolevels['cells'];
        sendData['messenger_obolevel'] = this.props.obolevels['messenger'];


        console.log("NI Sankey reload graph data");
        console.log(sendData);

        axios.post(config.getRestAddress() + "/sankey",sendData, config.axiosConfig)
        .then(function (response) {

            console.log("NI Sankey Chart data fetch")
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
            <D3SankeyChart graph={this.state.graphData} graphtitle={this.props.graphtitle}/>
        </div>);
    }

}