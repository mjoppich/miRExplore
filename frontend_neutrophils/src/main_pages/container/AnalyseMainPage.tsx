import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';

import QueryComponent from '../components/QueryComponent';

export interface AnalyseMainPageProps { };
export interface AnalyseMainPageState { queriesStored: number};

export class AnalyseMainPage extends React.Component<AnalyseMainPageProps, AnalyseMainPageState> {

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
        this.allQueries.push(
            <QueryComponent key={this.allQueries.length} loadSentences={false} showEvidenceTable={false}/>
        );
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

                        <p>In the analyse view you can query similarly to the explore View. However, results are shown in an aggregated manner.</p>
                        <p>You can create (interactive) plots, analyse the interaction network. But you will not be able to explore every result on an per-evidence base.</p>
                        <p>If you want to explore the network on a per-evidence base, you must switch into the <i>explore</i> mode.</p>

                    </CardText>
                </Card>
                
                <div>
                {this.allQueries}
                </div>

            </div>

        );
    }

}
