import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';

import Paper from 'material-ui/Paper';

//import {cytoscape} from 'cytoscape';
var cytoscape = require('cytoscape');

let cyStyle = {
  height: '400px',
  display: 'block'
};

export interface CytoscapeProps { elements: any };
export interface CytoscapeState { queriesStored: number};

class Cytoscape extends React.Component<CytoscapeProps, CytoscapeState> {
  cy = null;

  componentDidMount(){

    let cy = cytoscape({style: cyStyle, container: this.refs.cyelement});

    this.cy = cy;
    cy.json({elements: this.props.elements});
  }

  shouldComponentUpdate(){
    return false;
  }

  componentWillReceiveProps(nextProps){
    this.cy.json(nextProps);
  }

  componentWillUnmount(){
    this.cy.destroy();
  }

  getCy(){
    return this.cy;
  }

  render(){
    return <div style={cyStyle} ref="cyelement" />
  }
}


export interface CyParentProps { };
export interface CyParentState { };
class CyParent extends React.Component<CyParentProps, CyParentState>{
    componentDidMount(){
      // this is a good place for events
      //this.refs.graph.getCy();
    }
  
    render(){
      return (
          <Cytoscape ref="graph" elements={[{data: { id: 'a' }}]} />
      )
    }
}
  
export interface QueryComponentProps { key: number};
export interface QueryComponentState { };
class QueryComponent extends React.Component<QueryComponentProps, QueryComponentState> {
    constructor(props) {
        super(props);

    }

    render()
    {
        return (

            <Card>
                <CardHeader
                title="Query"
                subtitle="Subtitle"
                />
                <CardText>

                    <p>This is a query!</p>
                    <CyParent/>

                </CardText>
            </Card>

        );
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
        console.log("Query added");
    }

    clearQueries()
    {
        this.allQueries = [];
        this.setState({queriesStored: this.allQueries.length});
        console.log("Queries cleared");
    }

    /**
     * Render the component.
     */
    render() {

        console.log("ExploreMainpage render");


        return (

            <Paper>

                <Card>
                    <CardHeader
                    title="Without Avatar"
                    subtitle="Subtitle"
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

            </Paper>

        );
    }

}