import * as React from 'react'; 

import ChipInput from 'material-ui-chip-input'
import AutoComplete from 'material-ui/AutoComplete'
import axios from 'axios';
import config from '../config';

import FlatButton from 'material-ui/FlatButton';
import CheckIcon from 'material-ui/svg-icons/action/done';
import DisagreeIcon from 'material-ui/svg-icons/action/thumb-down';


interface IEntityElement {
    name: string;
    group: string;
 }

export interface EvidenceReportButtonProps { onAccept: any, onDisagree: any, dataID: any};
export interface EvidenceReportButtonState { acceptBgColor:string, disagreeBgColor: string, fdata: any};


export default class EvidenceReportButtons extends React.Component<EvidenceReportButtonProps, EvidenceReportButtonState>{

    constructor(props)
    {
        super(props);
    }

    componentWillMount()
    {
        var self=this;

        this.setState({acceptBgColor: "", disagreeBgColor: "", fdata: null});
    }

    componentDidMount(){
        this.updateColors();
    }

    updateColors()
    {

        var self = this;


        if (this.props.dataID)
        {
            var sendData = {'data_id': this.props.dataID}

            axios.post(config.getRestAddress() + "/feedback_info",sendData, config.axiosConfig)
            .then(function (response) {
              console.log(response.data);
    
              if (response.data.feedback_cum > 0)
              {
                  self.setState({acceptBgColor: "#D0E2BF", disagreeBgColor: "", fdata: response.data});

              } else  if (response.data.feedback_cum < 0)
              {
                self.setState({acceptBgColor: "", disagreeBgColor: "#e2bfbf", fdata: response.data});

              } else
              {
                self.setState({acceptBgColor: "", disagreeBgColor: "", fdata: response.data});
              }
    
            })
            .catch(function (error) {
              console.log(error);
            });
        }
    }

    onFeedback(accept)
    {

        if (accept)
        {
            this.props.onAccept();
        } else {
            this.props.onDisagree();
        }

        this.updateColors();
    }
  
    render(){


        return (
            <div>
                <FlatButton label="Accept" backgroundColor={this.state.acceptBgColor} onClick={() => this.onFeedback(true)} icon={<CheckIcon/>}/>
                <FlatButton label="Disagree" backgroundColor={this.state.disagreeBgColor} onClick={() => this.onFeedback(false)} icon={<DisagreeIcon/>}/>
            </div>
            );

    }
}