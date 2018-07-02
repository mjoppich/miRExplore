import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';

import QueryComponent from '../components/QueryComponent';
import FlatOboSynViewer from '../components/FlatOboSynViewer';

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


                <Card style={{marginBottom: "20px"}}>
                    <CardHeader
                    title="Cells"
                    subtitle="Lists all cells"
                    />

                    <CardText>

                        <p>This ontology currently reflects the <a target="_blank" href="https://www.ebi.ac.uk/ols/ontologies/cl">Cell Ontology</a>.</p>
                        <p>Currently, hierarchies are first generated from the <a target="_blank" href="https://www.ebi.ac.uk/ols/ontologies/cl/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FCL_0000255">
eukaryotic cell</a> obo term, then from any remaining term.</p>

                        <FlatOboSynViewer oboname="cells"/>

                    </CardText>
                </Card>

                <Card style={{marginBottom: "20px"}}>
                    <CardHeader
                    title="Messengers"
                    subtitle="Lists all messengers"
                    />

                    <CardText>

                        <p>This list has been automatically generated and needs manual curation. This is a flat hierarchy which is generally not preferable.</p>

                        <FlatOboSynViewer oboname="messengers"/>

                    </CardText>
                </Card>

                <Card style={{marginBottom: "20px"}}>
                    <CardHeader
                    title="Categories"
                    subtitle="Lists all categories"
                    />

                    <CardText>
                        <p>This list has been manually generated. This is a flat hierarchy which is generally not preferable.</p>
                        <FlatOboSynViewer oboname="categories"/>
                    </CardText>
                </Card>

            </div>

        );
    }

}
