import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';


export interface AnalyseMainPageProps { switchTab?: any };
export interface AnalyseMainPageState { };

export class AnalyseMainPage extends React.Component<AnalyseMainPageProps, AnalyseMainPageState> {

    allQueries = [];
    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);
    }


    /**
     * Render the component.
     */
    render() {


        return (

            <div>

                <Card style={{marginBottom: "20px"}}>
                    <CardHeader
                    title="Input Data Description"
                    subtitle="This is with what we work"
                    />
                    <CardActions>

                    </CardActions>
                    <CardText>

                        <p>For any bioinformatics analysis it is essential to know from which data results have been derived.</p>
                        <p>Here the used ontologies and the included synonym/search words are presented. For this analysis, Pubmed abstracts (last upd. January 2018) have been searched for the below listed terms.</p>
                        <p>As always, this analysis is not static. Many ideas can be implemented. When the results are almost correct, updating the pubmed database is anticipated.</p>
                        <p>If you want to explore the network on a per-evidence base, you must switch into the <i>explore</i> mode (bottom).</p>

                    </CardText>
                </Card>


            </div>

        );
    }

}
