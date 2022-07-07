import * as React from 'react'; 

import axios from 'axios';
import config from '../config';

import Button from '@mui/material/Button'
import CheckCircle from '@mui/icons-material/CheckCircle'
import Cancel from '@mui/icons-material/Cancel'



export interface EvidenceReportButtonProps { onAccept: any, onDisagree: any, dataID: any};
export interface EvidenceReportButtonState { acceptBgColor:string, disagreeBgColor: string, fdata: any};


export default class EvidenceReportButtons extends React.Component<EvidenceReportButtonProps, EvidenceReportButtonState>{

    constructor(props:any)
    {
        super(props);
        this.state = {acceptBgColor: "secondary", disagreeBgColor: "secondary", fdata: null};
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
    
              if (response.data.feedback_cum > 0)
              {
                  self.setState({acceptBgColor: "success", disagreeBgColor: "secondary", fdata: response.data});

              } else  if (response.data.feedback_cum < 0)
              {
                self.setState({acceptBgColor: "secondary", disagreeBgColor: "error", fdata: response.data});

              } else
              {
                self.setState({acceptBgColor: "secondary", disagreeBgColor: "secondary", fdata: response.data});
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

        var self=this;

        return (
            <div>
                <Button color={self.state.acceptBgColor} onClick={() => self.onFeedback(true)} startIcon={<CheckCircle/>}>Support</Button>
                <Button color={self.state.disagreeBgColor} onClick={() => self.onFeedback(false)} startIcon={<Cancel/>}>Disagree</Button>
            </div>
            );

    }
}