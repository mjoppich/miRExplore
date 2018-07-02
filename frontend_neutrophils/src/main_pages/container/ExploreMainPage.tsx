import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';

import QueryComponent from '../components/QueryComponent';

export interface ExplorePageProps { switchTab?: any };
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
                    subtitle="Explore new interactions"
                    />
                    <CardActions>
                    <FlatButton label="New Query" onClick={this.newQuery.bind(this)}/>
                    <FlatButton label="Clear" onClick={this.clearQueries.bind(this)}/>
                    </CardActions>
                    <CardText>

                        <p>Using a query you can query our database for interactions.</p>
                        <p>You must select an entity for which you would like to explore interactions. This might be one term, but can also be multiple terms.</p>
                        <p>You can restrict your search to such evidences, which must be for a specific organism, or must include a specific messenger or category.</p>
                        <p>You can show a interaction network for your selected evidences, or a sankey chart, showing you also how the interactions are connected to messengers or categories</p>
                        <p>Enjoy exploring your neutrophils!</p>


                    </CardText>
                </Card>
                
                <div>
                {this.allQueries}
                </div>

            </div>

        );
    }

}
