import * as React from "react"; 
import * as ReactDOM from "react-dom";
import FeatureViewer from './FeatureViewer';

import * as pyc from "pycollections";


export interface MEFeatureViewerProps { features: any, location:any}
export interface MEFeatureViewerState { }



export default class MEFeatureViewer extends React.Component<MEFeatureViewerProps, MEFeatureViewerState> {

    inputJSON: any = null;


    constructor(props) {
        super(props);

        console.log(this.props);

        this.inputJSON =  [{
            "gene_id": "NONMMUG000001.2",
            "gene_type": "lncrna",
            "gene_length": 2000,
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
                    "mirnas": {"mir103": [{"x": 15, "y": 35}]}
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
        }];


    }

    randomSeq(seqLen: number) {
        var text = "";
        var possible = "ABCDEFGHIJKLMNOPQRSTUVWXY";
      
        for (var i = 0; i < seqLen; i++)
          text += possible.charAt(Math.floor(Math.random() * possible.length));
      
        return text;
    }
      

    componentDidMount()
    {
        const de = ReactDOM.findDOMNode(this.refs.fv1);

        const jqelem = $(de);

        console.log("fv mount")
        console.log(de);
        
        for (var i = 0; i < this.inputJSON.length; ++i)
        {
            this.makeFeatureViewer(de, this.inputJSON[i]);
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

         var ft = new FeatureViewer(this.randomSeq(geneLength), component, {
            showAxis: true,
            showSequence: true,
            brushActive: true, //zoom
            toolbar:true, //current zoom & mouse position
            bubbleHelp:true, 
            zoomMax:50 //define the maximum range of the zoom
        });

         var transcripts = featureJSON['transcripts_list'];
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
                color: "#378a72",
                type: "rect",
                filter: "type1"
            });
         }



    }

    render() {

        console.log("ME FeatureViewer Render");

        return (<div><div ref="fv1" style={{width: "1000px"}}></div></div>);
    }

}