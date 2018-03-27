import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';



export class MainStatus extends React.Component<{},{}> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

    }

    /**
     * This method will be executed after initial rendering.
     */
    componentDidMount() {

    }

    /**
     * Render the component.
     */
    render() {


        return (
            <div>
                <p>
                    <h4>General Information</h4>
                    <ul>
                        <li>Number of Genes</li>
                        <li>Number of miRNA</li>
                        <li>Number of Interactions</li>
                    </ul>
                </p>
                <p>
                    <h4>Included Databases</h4>
                    <ul>
                        <li>miRWalk</li>
                        <li>miRecords</li>
                        <li>miDIP</li>
                    </ul>
                </p>
            </div>
        );
    }

}


export class WelcomePage extends React.Component<{},{}> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

    }

    /**
     * This method will be executed after initial rendering.
     */
    componentDidMount() {

    }

    /**
     * Render the component.
     */
    render() {


        return (

            <div>
            <Card style={{marginBottom: "20px"}}>
                <CardTitle
                    title="Welcome!"
                    subtitle="miRExplore v0.1"
                />
                <CardText >
                    <p>miRExplore v.01 online.</p>
                    <p>Explore the database in the explore tab</p>
                    <p>Provide some statistics of the project here</p>

                </CardText>
            </Card>
                        <Card>
                        <CardTitle
                            title="miRExplore"
                            subtitle="Statistics"
                        />
                        <CardText >
                            <MainStatus/>
                        </CardText>
                    </Card>
            </div>

        );
    }

}