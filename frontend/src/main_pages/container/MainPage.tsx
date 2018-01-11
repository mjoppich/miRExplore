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
                    title=""
                />
                <CardText >
                <img style={{width: 480+'px', display: "block", margin: "auto"}} src={`/img/welcome.png`} />
                </CardText>
            </Card>

        );
    }

}