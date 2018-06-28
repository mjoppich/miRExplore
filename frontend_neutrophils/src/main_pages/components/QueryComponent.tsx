import * as React from 'react';
import Switch from '@material-ui/core/Switch'
import OrganismChipAC from '../components/OrganismChipAC';
import EntityChipAC from '../components/EntityChipAC';
import LinearProgress from '@material-ui/core/LinearProgress';
import OboChipAC from '../components/OBOChipAC';
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import Button from '@material-ui/core/Button';
import axios from 'axios';
import config from '../config';
import TextField from '@material-ui/core/TextField';

import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel'

import QueryResult from './QueryResult';

export interface QueryComponentProps { key: number, loadSentences?: boolean, showEvidenceTable?: boolean, showShortEvidenceTable?: boolean};
export interface QueryComponentState {
    selectedElements: Array<any>,
    selectedOrganisms: Array<any>,
    selectedCategories: Array<any>,
    selectedMessengers: Array<any>,
    interactions: any,

    showInteractionGraph:boolean,
    showSankeyChart: boolean,
    queryStarted: boolean,
    obolevel: number,
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
        selectedElements: [],
        selectedOrganisms: [],
        selectedMessengers: [],
        selectedCategories: [],
        interactions: [],

        showSankeyChart:true,
        showInteractionGraph: false,

        obolevel: 3,
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
            return;
        }

        for (var i = 0; i < this.state.selectedElements.length; ++i)
        {
            sendData['elements'] = this.state.selectedElements;
        }

        if ((this.state.selectedMessengers) && (this.state.selectedMessengers.length > 0))
        {
            sendData['messengers'] = this.state.selectedMessengers;
        }

        if ((this.state.selectedCategories) && (this.state.selectedCategories.length > 0))
        {
            sendData['categories'] = this.state.selectedCategories;
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
                alignResults.push(<p key={alignResults.length}>No Result Available for your query.</p>)   ;
            }
            //alignResults.push(<pre key={1}>{JSON.stringify(this.state.selectedOrganisms, null, 2)}</pre>)   ;

        } else {

            alignResults.push(
                <QueryResult
                key={alignResults.length}
                showEvidenceTable={this.props.showEvidenceTable}
                showShortEvidenceTable={this.props.showShortEvidenceTable}
                showInteractionGraph={this.state.showInteractionGraph}
                showSankeyChart={this.state.showSankeyChart}
                searchWords={this.state.selectedElements}
                foundRelations={this.state.interactions["rels"]}
                docInfos={this.state.interactions["pmidinfo"]}
                searchQuery={this.state.interactions["searchquery"]}
                obolevel={this.state.obolevel}
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
                            floatText="Cell-Type"
                            hintText="Enter Cell-Type"
                            onValueChange={(newvalues) => this.setState({selectedElements: newvalues})
                        }/>

                        <OboChipAC
                            url="organisms"
                            floatText="Organism"
                            hintText="Enter organism name"
                            onValueChange={(newvalues) => this.setState({selectedOrganisms: newvalues})
                        }/>


                        <OboChipAC
                            url="category_ac"
                            floatText="Categories"
                            hintText="Enter category name"
                            onValueChange={(newvalues) => this.setState({selectedCategories: newvalues})
                        }/>

                        <OboChipAC
                            url="messengers_ac"
                            floatText="Messengers"
                            hintText="Enter Messenger-Term here"
                            onValueChange={(newvalues) => this.setState({selectedMessengers: newvalues})
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
                checked={this.state.showSankeyChart}
                onChange={(newValue) => this.setState({showSankeyChart: !this.state.showSankeyChart})}
                /> 
            }
            label="Show Sankey Chart"
            />

            <TextField
                id="oboinput"
                label="Ontology Hierarchy Level"
                value={this.state.obolevel}
                onChange={(evt) => {console.log(evt.target.value); this.setState({obolevel: Number(evt.target.value)})}}
                type="number"
                helperText="For plots only: specify how many levels to go down in cell hierarchy"

                InputLabelProps={{
                shrink: true,

                }}
                margin="normal"
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
