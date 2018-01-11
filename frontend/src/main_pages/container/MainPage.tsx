import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';

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

            <Card>
                <CardTitle
                    title="Welcome!"
                    subtitle="miRExplore v0.1"
                />
                <CardText >
                    <p>miRExplore v.01 now online.</p>
                    <p>Following PubMed state: 2018 base</p>

                </CardText>
            </Card>

        );
    }

}