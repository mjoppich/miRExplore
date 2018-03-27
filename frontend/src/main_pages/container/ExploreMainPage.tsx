import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';
import Paper from 'material-ui/Paper';
import SelectedElements from '../components/SelectedElements';
import ACInput from '../components/AutoComplete';
import MSATableViewer from '../components/MSATableViewer';
import axios from 'axios';
import config from '../config';
import D3GraphViewer from '../components/D3GraphViewer';
import Toggle from 'material-ui/Toggle';
import OrganismChipAC from '../components/OrganismChipAC';
import EntityChipAC from '../components/EntityChipAC';

import * as Collections from 'typescript-collections';

export interface D3ParentProps { graphData: any};
export interface D3ParentState { graph: any };
class D3GraphComponent extends React.Component<D3ParentProps, D3ParentState>{

    neo4jd3: any = null;

    constructor(props)
    {
        super(props);

    }

    componentWillMount()
    {
    }

    componentDidMount(){

    }
 
    render(){
      return (
          <D3GraphViewer id="bla" />
      )
    }
}


  
export interface QueryComponentProps { key: number};
export interface QueryComponentState { 
    selectedElements: Array<any>,
    selectedOrganisms: Array<any>,
    interactions: any,
    matureMIRNA: boolean,
    restrictOrganism: boolean
 };

class QueryComponent extends React.Component<QueryComponentProps, QueryComponentState> {
    constructor(props) {
        super(props);

    }

    componentWillMount()
    {
        this.setState({selectedElements: [], matureMIRNA:true, restrictOrganism: false});
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

        console.log("Selected Elements in EXplore")
        console.log(this.state.selectedElements);
        console.log(sendData)

        for (var i = 0; i < this.state.selectedElements.length; ++i)
        {
            let elem = this.state.selectedElements[i];

            let elemGroup = elem['group'];
            let elemName = elem['name'];

            console.log(elemGroup)
            console.log(elemName)

            if (elemGroup in sendData)
            {
                sendData[elemGroup].push(elemName);
            } else {
                sendData[elemGroup] = [elemName];
            }
        }

        console.log(sendData);

        axios.post(config.getRestAddress() + "/find_interactions",sendData, config.axiosConfig)
        .then(function (response) {
          console.log(response.data)

          self.setState({interactions: response.data})

        })
        .catch(function (error) {
          console.log(error);
          self.setState({interactions: {}});
        });
    }



    render()
    {

        var alignResults = [];
        if ((this.state.interactions == null) || (this.state.interactions.length == 0))
        {
            alignResults.push(<p key={0}>No Result Available yet</p>)   ;
            alignResults.push(<pre key={1}>{JSON.stringify(this.state.selectedOrganisms, null, 2)}</pre>)   ;

        } else {

            alignResults.push(<pre key={0}>{JSON.stringify(this.state.interactions, null, 2)}</pre>)   ;
            alignResults.push(<pre key={1}>{JSON.stringify(this.state.selectedOrganisms, null, 2)}</pre>)   ;

            //var alignKeys = Object.keys(this.state.alignments);        
        }

        return (<Card style={{marginBottom: "20px"}}>
                <CardHeader
                title="Search Homology Entries"
                subtitle="Search by Gene/Protein ID"
                />
                <CardText>

                    <div>
                        <EntityChipAC onValueChange={(newvalues) => {console.log("onVC called"); this.setState({selectedElements: newvalues})}} />
                        <OrganismChipAC onValueChange={(newvalues) => this.setState({selectedOrganisms: newvalues})}/>


                <Toggle
                label="Use mature miRNA (instead of pre-miRNA)"
                defaultToggled={true}
                toggled={this.state.matureMIRNA}
                onToggle={(event, newValue) => this.setState({matureMIRNA: !this.state.matureMIRNA})}
                />
                <Toggle
                label="Restrict to specified organisms"
                defaultToggled={false}
                toggled={this.state.restrictOrganism}
                onToggle={(event, newValue) => this.setState({restrictOrganism: !this.state.restrictOrganism})}
                />
                        <FlatButton label="Query specified Elements" onClick={() => this.prepareResults()}/>
                    </div>

                    <div>
                        {alignResults}
                    </div>

                </CardText>
            </Card>);
    }
};

export interface ExplorePageProps { };
export interface ExplorePageState { queriesStored: number};

export class ExploreMainPage extends React.Component<ExplorePageProps, ExplorePageState> {

    allQueries = [];
    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);
    }

    /**
     * This method will be executed after initial rendering.
     */
    componentWillMount() {

        this.setState({queriesStored: 0});

    }

    newQuery()
    {
        this.allQueries.push(<QueryComponent key={this.allQueries.length}/>);
        this.setState({queriesStored: this.allQueries.length});
    }

    clearQueries()
    {
        this.allQueries = [];
        this.setState({queriesStored: this.allQueries.length});
    }

    /**
     * Render the component.
     */
    render() {


        return (

            <div>

                <Card style={{marginBottom: "20px"}}>
                    <CardHeader
                    title="Create Query"
                    subtitle="Manage your queries"
                    actAsExpander={true}
                    showExpandableButton={true}
                    />
                    <CardActions>
                    <FlatButton label="New Query" onClick={this.newQuery.bind(this)}/>
                    <FlatButton label="Clear" onClick={this.clearQueries.bind(this)}/>
                    </CardActions>
                    <CardText expandable={true}>

                        <p>Some explanation on how to use this view!</p>

                    </CardText>
                </Card>
                
                <div>
                {this.allQueries}
                </div>

            </div>

        );
    }

}