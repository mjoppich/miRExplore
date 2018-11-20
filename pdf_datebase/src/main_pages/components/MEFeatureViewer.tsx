import * as React from "react"; 
import * as ReactDOM from "react-dom";
import FeatureViewer from './FeatureViewer';
import axios from 'axios';
import config from '../config';

import * as pyc from "pycollections";


export interface MEFeatureViewerProps { featureID: string}
export interface MEFeatureViewerState { features: any, neighbourhood:any }

export default class MEFeatureViewer extends React.Component<MEFeatureViewerProps, MEFeatureViewerState> {

    fv1: any = null;
    fv2: any = null;

    readonly state = {
        'features': [],
        'neighbourhood': {}
    };


    constructor(props) {
        super(props);

        this.fv1 = React.createRef();
        this.fv2 = React.createRef();

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

    getNewNeighbourhood()
    {
        console.log("Attempting to fetch neighbours");


        var self=this;
        axios.get(config.getRestAddress() + "/get_neighbours/"+this.props.featureID, config.axiosConfig)
        .then(function (response) {

            console.log("Feature Viewer Data Receive: neighbours")
            console.log(response.data)

            self.setState({neighbourhood: response.data});
            
        })
        .catch(function (error) {
          console.log(error);
          self.setState({neighbourhood: {}});
        });
    }

    componentDidUpdate(prevProps, prevState, snapshot) {

        var self=this;

        if (prevProps != this.props)
        {

            this.getNewFeatures();
            this.getNewNeighbourhood();

        }

        if (this.state.features != prevState.features)
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

        if (this.state.neighbourhood != prevState.neighbourhood)
        {
            var domNode = this.fv2.current as HTMLElement;

            while (domNode.firstChild)
            {
                domNode.removeChild(domNode.firstChild);
            }

            var allOrgs = Object.keys(this.state.neighbourhood);

            for (var i = 0; i < allOrgs.length; ++i)
            {


                var allgenes = Object.keys(this.state.neighbourhood[allOrgs[i]]);

                for (var j = 0; j < allgenes.length; ++j)
                {
                    var newchild = document.createElement("div");
                    domNode.appendChild(newchild);


    
                    this.makeNeighbourViewer(newchild, this.state.neighbourhood[allOrgs[i]][allgenes[j]]);
                }
            }
        }

    }

    makeNeighbourViewer(component, featureJSON)
    {
        /**
         * 
         * 
         *
         
{
  "ENSMUSG00000025907": {
    "ENSMUSG00000025907": {
      "chr": "chr1",
      "start": 6206197,
      "end": 6276648,
      "strand": "+",
      "show_start": 42124,
      "show_end": 112575
    },
    "nb": [
      {
        "chr": "chr1",
        "start": 6204110,
        "end": 6205521,
        "strand": "-",
        "name": "LNC_GE_mm10_00536088",
        "show_start": 40037,
        "show_end": 41448
      },
      {
        "chr": "chr1",
        "start": 6183215,
        "end": 6184262,
        "strand": "-",
        "name": "LNC_GE_mm10_00536081",
        "show_start": 19142,
        "show_end": 20189
      },
      {
        "chr": "chr1",
        "start": 6180675,
        "end": 6182236,
        "strand": "-",
        "name": "LNC_GE_mm10_00536080",
        "show_start": 16602,
        "show_end": 18163
      },
      {
        "chr": "chr1",
        "start": 6288197,
        "end": 6289173,
        "strand": "+",
        "name": "LNC_GE_mm10_00534657",
        "show_start": 124124,
        "show_end": 125100
      },
      {
        "chr": "chr1",
        "start": 6209866,
        "end": 6215293,
        "strand": "-",
        "name": "LNC_GE_mm10_00496805",
        "show_start": 45793,
        "show_end": 51220
      },
      {
        "chr": "chr1",
        "start": 6191675,
        "end": 6192489,
        "strand": "-",
        "name": "LNC_GE_mm10_00536083",
        "show_start": 27602,
        "show_end": 28416
      },
      {
        "chr": "chr1",
        "start": 6283373,
        "end": 6287775,
        "strand": "+",
        "name": "LNC_GE_mm10_00534656",
        "show_start": 119300,
        "show_end": 123702
      },
      {
        "chr": "chr1",
        "start": 6192603,
        "end": 6193379,
        "strand": "-",
        "name": "LNC_GE_mm10_00536084",
        "show_start": 28530,
        "show_end": 29306
      },
      {
        "chr": "chr1",
        "start": 6176944,
        "end": 6178087,
        "strand": "-",
        "name": "LNC_GE_mm10_00536079",
        "show_start": 12871,
        "show_end": 14014
      },
      {
        "chr": "chr1",
        "start": 6194082,
        "end": 6195121,
        "strand": "-",
        "name": "LNC_GE_mm10_00536085",
        "show_start": 30009,
        "show_end": 31048
      },
      {
        "chr": "chr1",
        "start": 6189484,
        "end": 6190754,
        "strand": "-",
        "name": "LNC_GE_mm10_00536082",
        "show_start": 25411,
        "show_end": 26681
      },
      {
        "chr": "chr1",
        "start": 6195273,
        "end": 6196193,
        "strand": "-",
        "name": "LNC_GE_mm10_00536086",
        "show_start": 31200,
        "show_end": 32120
      },
      {
        "chr": "chr1",
        "start": 6196344,
        "end": 6203804,
        "strand": "-",
        "name": "LNC_GE_mm10_00536087",
        "show_start": 32271,
        "show_end": 39731
      },
      {
        "chr": "chr1",
        "start": 6206197,
        "end": 6276648,
        "strand": "+",
        "name": "ENSMUSG00000025907",
        "show_start": 42124,
        "show_end": 112575
      },
      {
        "chr": "chr1",
        "start": 6274564,
        "end": 6276365,
        "strand": "+",
        "name": "LNC_GE_mm10_00496806",
        "show_start": 110491,
        "show_end": 112292
      },
      {
        "chr": "chr1",
        "start": 6287667,
        "end": 6288962,
        "strand": "-",
        "name": "LNC_GE_mm10_00536089",
        "show_start": 123594,
        "show_end": 124889
      },
      {
        "chr": "chr1",
        "start": 6164173,
        "end": 6167111,
        "strand": "-",
        "name": "LNC_GE_mm10_00496804",
        "show_start": 100,
        "show_end": 3038
      },
      {
        "chr": "chr1",
        "start": 6173972,
        "end": 6176068,
        "strand": "-",
        "name": "LNC_GE_mm10_00536078",
        "show_start": 9899,
        "show_end": 11995
      },
      {
        "chr": "chr1",
        "start": 6359218,
        "end": 6394731,
        "strand": "+",
        "name": "ENSMUSG00000087247",
        "show_start": 195145,
        "show_end": 230658
      }
    ],
    "range": {
      "start": 0,
      "end": 230658
    }
  }
}



         * 
         * 
         */

         console.log("FeatureJSON")
         console.log(featureJSON)

        var targetGenes = Object.keys(featureJSON);

        var removeKeys = ['nb', 'range'];

        var rangeStart = Number(featureJSON['range']['start']);
        var rangeEnd = Number(featureJSON['range']['end']);



        var ft = new FeatureViewer(rangeEnd-rangeStart, component, {
            showAxis: true,
            showSequence: true,
            sequenceName: "",
            brushActive: true, //zoom
            toolbar:true, //current zoom & mouse position
            bubbleHelp:true, 
            zoomMax:50 //define the maximum range of the zoom
        });

        var compareByName = function(a,b) {
            if (a.name < b.name)
              return -1;
            if (a.name > b.name)
              return 1;
            return 0;
        }

        targetGenes = targetGenes.sort(compareByName)


        for (var i=0; i < targetGenes.length; ++i)
        {

            if (removeKeys.indexOf(targetGenes[i]) >= 0)
            {
                console.log("Ignoring targetGenes " + targetGenes[i])
                continue;
            }

            var tgElem = featureJSON[targetGenes[i]];

            ft.addFeature({
                data: [{
                    x: tgElem['show_start'],
                    y: tgElem['show_end'],
                    id: tgElem['name'],
                    description: tgElem['name']
                }],
                name: tgElem['name'],
                className: "test1",
                color: "#378a72",
                type: "rect",
                filter: "type1"
            });

        }

        var allNeighbours =  featureJSON['nb'];

        compareByName = allNeighbours.sort(compareByName);

        for (var i=0; i <allNeighbours.length; ++i)
        {

            var tgElem = allNeighbours[i];

            ft.addFeature({
                data: [{
                    x: tgElem['show_start'],
                    y: tgElem['show_end'],
                    id: tgElem['name'],
                    description: tgElem['name']
                }],
                name: tgElem['name'],
                className: "test1",
                color: "#779000",
                type: "rect",
                filter: "type1"
            });

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

         var allMirnaBS = new Array<any>();
         var addedMirnaNames = new Array<any>();

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

                console.log("Total amount of mirnas:");
                console.log(relKeys.length)

                var allMirnaLengths = relKeys.length;

                for (var m= 0; m < relKeys.length; ++m)
                {
                    var mirnaID = relKeys[m];
                    var mirnaBS = allMirnas[ mirnaID ];

                    if (addedMirnaNames.indexOf(mirnaID) >= 0)
                    {
                        continue;
                    }

                    // do these mirnaBS fit into a single existing bin?
                    var added = false;

                    if (allMirnaLengths > 30)
                    {
                        for (var ibins = 0; ibins < allMirnaBS.length; ++ibins)
                        {
                            var fitInBin = true;
                            var bin = allMirnaBS[i];
    
                            for (var ibs=0; ibs < mirnaBS.length; ++ibs)
                            {
                                
                                for (var ibinElem; ibinElem < bin.length; ++ibinElem)
                                {
                                    if ( (( elem.x <= mirnaBS[ibs].x) && (mirnaBS[ibs].x <= elem.y)) || (( elem.x <= mirnaBS[ibs].y) && (mirnaBS[ibs].y <= elem.y  )))
                                    {
                                        fitInBin = false;
                                    }
                                }
                            }
    
                            if (fitInBin)
                            {
                                for (var ibs = 0; ibs < mirnaBS.length; ++ibs)
                                {
                                    var element = mirnaBS[i];
                                    bin.push( 
                                        {
                                            x: element.x,
                                            y: element.y,
                                            id: mirnaID,
                                            description: mirnaID
                                        }
                                     )
                                }
                                added=true;
    
                                addedMirnaNames.push(mirnaID);
    
                                allMirnaBS[i] = bin;
    
                                break;
                            }
    
    
                        }
                    }
                    

                    if (!added)
                    {
                        var newbin = new Array<any>();
                        for (var ibs = 0; ibs < mirnaBS.length; ++ibs)
                        {
                            var element = mirnaBS[i];
                            newbin.push( 
                                {
                                    x: element.x,
                                    y: element.y,
                                    id: mirnaID,
                                    description: mirnaID
                                }
                             )
                        }
                        addedMirnaNames.push(mirnaID);

                        allMirnaBS.push(newbin);
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

        console.log("Need to add bins")
        console.log(allMirnaBS.length);

         for (var i = 0; i < allMirnaBS.length; ++i)
         {

            var allNames = [];
            var allPositions = []

            for (var j = 0; j < allMirnaBS[i].length; ++j)
            {
                var elem = allMirnaBS[i][j];
                var elemID = elem.id;

                if (allNames.indexOf(elemID) < 0)
                {
                    allNames.push(elemID);
                }

                allPositions.push(elem);
                
            }

            ft.addFeature({
                data: allPositions,
                name: allNames.join(", "),
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
        name: "Rfam",
        className: "test1",
        color: "#FF8a72",
        type: "rect",
        filter: "type1"  
    })


    }

    componentDidMount()
    {
        this.getNewFeatures();
        this.getNewNeighbourhood();
    }

    render() {

        console.log("ME FeatureViewer Render");

        return (<div>
                    <div ref={this.fv1} style={{width: "1000px"}}></div>
                    <div ref={this.fv2} style={{width: "1000px"}}></div>                    
                </div>);
    }

}