import * as React from 'react';
import Switch from '@material-ui/core/Switch'
import LinearProgress from '@material-ui/core/LinearProgress';
import OboChipAC from '../components/OBOChipAC';
import {Card, CardHeader, CardText} from 'material-ui/Card';
import Button from '@material-ui/core/Button';
import axios from 'axios';
import config from '../config';

import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel'

import QueryResult from './QueryResult';

export interface QueryComponentProps { key: number, loadSentences?: boolean, showEvidenceTable?: boolean, showShortEvidenceTable?: boolean};
export interface QueryComponentState {
    selectedElements: Array<any>,
    selectedOrganisms: Array<any>,
    selectedDiseases: Array<any>,
    selectedGOs: Array<any>,
    selectedCells: Array<any>,
    selectedNCITs: Array<any>,

    interactions: any,

    showInteractionGraph:boolean,
    showFeatures: boolean,
    queryStarted: boolean,

    loadSentences: boolean
 };

export default class QueryComponent extends React.Component<QueryComponentProps, QueryComponentState> {

    public static defaultProps: Partial<QueryComponentProps> = {
        loadSentences: true,
        showEvidenceTable: true,
        showShortEvidenceTable: false
    };

    constructor(props) {
        super(props);

    }

    readonly state = {
        selectedElements:[],
        selectedOrganisms: [],
        selectedDiseases: [],
        selectedGOs: [],
        selectedCells: [],
        selectedNCITs: [],
    
        interactions: [],
    
        showInteractionGraph: false,
        showFeatures: false,
        queryStarted: false,
        loadSentences: null
    };

    static getDerivedStateFromProps(props, state) {
        if (state.loadSentences == null) {
          return {
            loadSentences: props.loadSentences,
          }
        }

        return null;
      }

    newElementSelected( newElement )
    {
        this.state.selectedElements.push(newElement);
        this.setState({selectedElements: this.state.selectedElements});
    }

    elementClicked( elementText )
    {
        console.log("Clicked on: " + elementText)
    }

    deleteElement( elementText, i )
    {
        var idx = this.state.selectedElements.indexOf(elementText);

        if (idx >= 0)
        {
            this.state.selectedElements.splice(idx, 1);
        }

        this.setState({selectedElements: this.state.selectedElements})
    }

    prepareResults()
    {
        var self = this;

        let sendData = {};

        console.log("Selected Elements in Explore")
        console.log(this.state.selectedElements);
        console.log(this.state);
        console.log(sendData);

        if (!(this.state.selectedElements) || (this.state.selectedElements.length == 0))
        {
            console.log("Empty selection of genes/mirnas")
            sendData["gene"] = [];
        } else {
            for (var i = 0; i < this.state.selectedElements.length; ++i)
            {
                let elem = this.state.selectedElements[i];
    
                let elemGroup = elem['group'];
                let elemName = elem['name'];
    
                if (elemGroup in sendData)
                {
                    sendData[elemGroup].push(elemName);
                } else {
                    sendData[elemGroup] = [elemName];
                }
            }
        }

        if ((this.state.selectedCells) && (this.state.selectedCells.length > 0))
        {
            sendData['cells'] = this.state.selectedCells;
        }

        if ((this.state.selectedDiseases) && (this.state.selectedDiseases.length > 0))
        {
            sendData['disease'] = this.state.selectedDiseases;
        }

        if ((this.state.selectedGOs) && (this.state.selectedGOs.length > 0))
        {
            sendData['go'] = this.state.selectedGOs;
        }

        if ((this.state.selectedNCITs) && (this.state.selectedNCITs.length > 0))
        {
            sendData['ncits'] = this.state.selectedNCITs;
        }

        if ((this.state.selectedOrganisms)&&(this.state.selectedOrganisms.length > 0))
        {
            sendData['organisms'] = this.state.selectedOrganisms;
        }

        if (!this.state.loadSentences)
        {
            sendData['sentences'] = "false";
        }

        console.log(sendData);

        this.setState({queryStarted: true});

        axios.post(config.getRestAddress() + "/find_interactions",sendData, config.axiosConfig)
        .then(function (response) {
          console.log(response.data)

          response.data.searchquery = sendData

          self.setState({interactions: response.data, queryStarted: false})

        })
        .catch(function (error) {
          console.log(error);
          self.setState({interactions: {}, queryStarted: false});
        });
    }



    render()
    {

/*


            <Toggle
            label="Include Predictive Interactions"
            defaultToggled={true}
            toggled={this.state.predictiveInteractions}
            onToggle={(event, newValue) => this.setState({predictiveInteractions: !this.state.predictiveInteractions})}
            />
            <Toggle
            label="Show Gene Structure"
            defaultToggled={false}
            toggled={this.state.showGeneStructure}
            onToggle={(event, newValue) => this.setState({showGeneStructure: !this.state.showGeneStructure})}
            />

                */


        var alignResults = [];

        if (this.state.queryStarted)
        {
            alignResults.push(<LinearProgress key={alignResults.length}/>);
            alignResults.push(<p key={alignResults.length}>Retrieving Results</p>);

        }

        if ((this.state.interactions == null) || (this.state.interactions.length == 0))
        {

            if (!this.state.queryStarted)
            {
                alignResults.push(<p key={alignResults.length}>Press Query Specified Elements to start exploring.</p>)   ;
            }
            //alignResults.push(<pre key={1}>{JSON.stringify(this.state.selectedOrganisms, null, 2)}</pre>)   ;

        } else {

            alignResults.push(
                <QueryResult
                key={alignResults.length}
                showEvidenceTable={this.props.showEvidenceTable}
                showShortEvidenceTable={this.props.showShortEvidenceTable}
                showInteractionGraph={this.state.showInteractionGraph}
                showFeatures={this.state.showFeatures}
                searchWords={this.state.selectedElements}
                foundRelations={this.state.interactions["rels"]}
                docInfos={this.state.interactions["pmidinfo"]}
                searchQuery={this.state.interactions["searchquery"]}
                
                />)

            //alignResults.push(<pre key={0}>{JSON.stringify(this.state.interactions, null, 2)}</pre>)   ;
            //alignResults.push(<pre key={1}>{JSON.stringify(this.state.selectedOrganisms, null, 2)}</pre>)   ;
            //var alignKeys = Object.keys(this.state.alignments);

            //<EntityChipAC onValueChange={(newvalues) => {console.log("onVC called"); this.setState({selectedElements: newvalues})}} />
        }

        return (<Card style={{marginBottom: "20px"}}>
                <CardHeader
                title="Search Interactions"
                subtitle="Browse by category/dimension"
                />
                <CardText>

                    <div>


                        <OboChipAC
                            url="autocomplete"
                            floatText="Gene Symbol/miRNA"
                            hintText="Enter Gene Symbol/miRNA"
                            onValueChange={(newvalues) => this.setState({selectedElements: newvalues})
                        }/>

                        <OboChipAC
                            url="organisms"
                            floatText="Organism"
                            hintText="Enter organism name"
                            onValueChange={(newvalues) => this.setState({selectedOrganisms: newvalues})
                        }/>

                        <OboChipAC
                            url="ncit_ac"
                            floatText="Protein Class"
                            hintText="Enter a Protein Class"
                            onValueChange={(newvalues) => this.setState({selectedNCITs: newvalues})
                        }/>

                        <OboChipAC
                            url="go_ac"
                            floatText="Gene Ontology"
                            hintText="Enter GO-Name"
                            onValueChange={(newvalues) => this.setState({selectedGOs: newvalues})
                        }/>

                        <OboChipAC
                            url="cells_ac"
                            floatText="Cells"
                            hintText="Enter Cell-name here"
                            onValueChange={(newvalues) => this.setState({selectedCells: newvalues})
                        }/>

                        <OboChipAC
                            url="disease_ac"
                            floatText="Diseases"
                            hintText="Enter Disease-name here"
                            onValueChange={(newvalues) => this.setState({selectedDiseases: newvalues})
                        }/>

        <FormGroup>

            
            <FormControlLabel
            control={
                <Switch
                checked={this.state.loadSentences}
                onChange={(newValue) => this.setState({loadSentences: !this.state.loadSentences})}
                /> 
            }
            label="Load Sentences"
            />

            <FormControlLabel
            control={
                <Switch
                checked={this.state.showInteractionGraph}
                onChange={(newValue) => this.setState({showInteractionGraph: !this.state.showInteractionGraph})}
                />
            }
            label="Show Interaction Graph"
            />

            <FormControlLabel
            control={
                <Switch
                checked={this.state.showFeatures}
                onChange={(newValue) => this.setState({showFeatures: !this.state.showFeatures})}
                /> 
            }
            label="Show Gene Features"
            />

                <Button onClick={() => this.prepareResults()}>
                Query specified Elements
                </Button>
            </FormGroup>



                    </div>

                    <div>
                        {
                            alignResults
                        }

                    </div>

                </CardText>
            </Card>);
    }
};
