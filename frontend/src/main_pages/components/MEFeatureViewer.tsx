import * as React from "react"; 
import * as ReactDOM from "react-dom";
import FeatureViewer from './FeatureViewer';
import axios from 'axios';
import config from '../config';

import * as pyc from "pycollections";


export interface MEFeatureViewerProps { featureID: string}
export interface MEFeatureViewerState { features: any }

export default class MEFeatureViewer extends React.Component<MEFeatureViewerProps, MEFeatureViewerState> {

    fv1: any = null;

    readonly state = {
        'features': [],
    };


    constructor(props) {
        super(props);

        this.fv1 = React.createRef();

    }

    randomSeq(seqLen: number) {
        var text = "";
        var possible = "ABCDEFGHIJKLMNOPQRSTUVWXY";
      
        for (var i = 0; i < seqLen; i++)
          text += possible.charAt(Math.floor(Math.random() * possible.length));
      
        return text;
    }

    getNewFeatures()
    {
        var self=this;
        axios.get(config.getRestAddress() + "/get_features/"+this.props.featureID, config.axiosConfig)
        .then(function (response) {

            console.log("Feature Viewer Data Receive")
            console.log(response.data)

            self.setState({features: response.data.features});
            
        })
        .catch(function (error) {
          console.log(error);
          self.setState({features: []});
        });
    }

    componentDidUpdate(prevProps, prevState, snapshot) {

        var self=this;

        if (prevProps != this.props)
        {

            this.getNewFeatures();

        }

        if (this.state.features)
        {
            var domNode = this.fv1.current as HTMLElement;

            while (domNode.firstChild)
            {
                domNode.removeChild(domNode.firstChild);
            }


            for (var i = 0; i < this.state.features.length; ++i)
            {
                var newchild = document.createElement("div");
                domNode.appendChild(newchild);

                console.log("Make Viewer")
                console.log(newchild);
                console.log(this.state.features[i]);

                this.makeFeatureViewer(newchild, this.state.features[i]);
            }
        }

    }
      
    makeFeatureViewer(component, featureJSON)
    {
        /**
         * 
         * 
         *
         {
            "gene_id": "NONMMUG000001.2",
            "gene_type": "lncrna",
            "transcripts_list": [
                {
                    "start": 0,
                    "stop": 100,
                    "coding_list": [
                        {
                            "x": 0,
                            "y": 20,
                            "type": "cds"
                        },
                        {
                            "x": 50,
                            "y": 70,
                            "type": "utr"
                        }            
                    ],
                    "transcript_id": "NONMMUT000001.2",
                    "mir103": [{"x": 15, "y": 35}]
                },
                {
                    "start": 0,
                    "stop": 200,
                    "coding_list": [
                        {
                            "x": 0,
                            "y": 50,
                            "type": "utr"
                        },
                        {
                            "x": 100,
                            "y": 130,
                            "type": "utr"
                        }            
                    ],
                    "transcript_id": "NONMMUT000001.1",
                    "mirnas": {"mir103": [{"x": 15, "y": 35}, {"x": 120, "y": 140}]}
                }       
            ]
        }
         * 
         * 
         */

         var geneID = featureJSON['gene_id'];
         var geneType = featureJSON['gene_type'];
         var geneLength = featureJSON['gene_length'];

         var geneStart = featureJSON['gene_start'];
         var geneStop = featureJSON['gene_stop'];
         var geneStrand = featureJSON['strand'];

        var rfams = featureJSON['rfams'];

         var geneSymbol = this.props.featureID;

         var transcripts = featureJSON['transcript_list'];

        if ((!transcripts) || (transcripts.length == 0))
        {
            return;
        }

         var ft = new FeatureViewer(geneLength, component, {
            showAxis: true,
            showSequence: true,
            sequenceName: geneID + "(" + geneSymbol +  ")",
            brushActive: true, //zoom
            toolbar:true, //current zoom & mouse position
            bubbleHelp:true, 
            zoomMax:50 //define the maximum range of the zoom
        });

         var allMirnaBS = {};

         for (var i = 0; i < transcripts.length; ++i)
         {
             var transcript = transcripts[i];

             /**
              *                 {
                    "start": 0,
                    "stop": 200,
                    "coding_list": [
                        {
                            "x": 0,
                            "y": 50,
                            "type": "utr"
                        },
                        {
                            "x": 100,
                            "y": 130,
                            "type": "utr"
                        }            
                    ],
                    "transcript_id": "NONMMUT000001.1",
                    "mirnas": {"mir103": [{"x": 15, "y": 35}, {"x": 120, "y": 140}]}
                } 
              */

                var codingData = transcript['coding_list'];
                codingData = codingData.map( (d) => {return {'x': d.x, 'y': d.y, 'description': d.type}; });

                var allMirnas = transcript['mirnas'];
                
                for (var attrname in transcript['other_mirnas'])
                {
                    allMirnas[attrname] = transcript['other_mirnas'][attrname];
                }

                console.log(codingData);
                console.log(transcript);

                var relKeys = Object.keys(allMirnas);

                for (var m= 0; m < relKeys.length; ++m)
                {
                    var mirnaID = relKeys[m];
                    var mirnaBS = allMirnas[ mirnaID ];

                    console.log("mirnaBS");
                    console.log(mirnaBS);

                    if (mirnaBS)
                    {
                        mirnaBS.forEach(element => {

                            if (mirnaID in allMirnaBS)
                            {
                                var idx = (allMirnaBS[mirnaID] as Array<any>).indexOf(element);
    
                                if (idx >= 0)
                                {
                                    allMirnaBS[mirnaID].push(element);
                                }
                            } else {
                                allMirnaBS[mirnaID] = [element];
                            }
    
    
                        });
                    }
                }

                var transcriptID = transcript['transcript_id'];

                ft.addFeature({
                    data: codingData,
                    name: transcriptID,
                    className: "test1",
                    color: "#378a72",
                    type: "rect",
                    filter: "type1"
                });


         }


         var selMirnaBS = Object.keys(allMirnaBS);

         for (var i = 0; i < selMirnaBS.length; ++i)
         {
            ft.addFeature({
                data: allMirnaBS[selMirnaBS[i]],
                name: selMirnaBS[i],
                className: "test1",
                color: "#378aFF",
                type: "rect",
                filter: "type1"
            });
         }


        var rfamFeatures = [];
         for (var i = 0; i < rfams.length; ++i)
         {
             // {'rfam_id': 'RF00005', 'orgid': 'mmu', 'chr': 'chr1', 'rfam_start': 4913784, 'rfam_end': 4913856, 'strand': '+', 'bit_score': 36.0, 'evalue_score': 0.44, 'cm_start': 1, 'cm_end': 71, 'truncated': 0, 'type': 'full'}
            var rfamElem = rfams[i];

            var rfamStart = rfamElem.rfam_start;
            var rfamEnd = rfamElem.rfam_end;
            var rfamName = rfamElem.rfam_id;

            if (geneStrand == '+')
            {
                rfamStart = rfamStart-geneStart;
                rfamEnd = rfamEnd-geneStart;
            } else { // neg strand
                var tmpStart = geneStop - rfamEnd;
                var tmpEnd = geneStop - rfamStart;

                rfamStart = tmpStart;
                rfamEnd = tmpEnd;
            }

            rfamFeatures.push({x: rfamStart, y: rfamEnd, description: rfamName});
         }


         ft.addFeature({
            data: rfamFeatures,
            name: selMirnaBS[i],
            className: "test1",
            color: "#FF8a72",
            type: "rect",
            filter: "type1"  
        })


    }

    componentDidMount()
    {
        this.getNewFeatures();
    }

    render() {

        console.log("ME FeatureViewer Render");

        return (<div><div ref={this.fv1} style={{width: "1000px"}}></div></div>);
    }

}