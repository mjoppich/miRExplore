


import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';

import FontIcon from 'material-ui/FontIcon';
import {BottomNavigation, BottomNavigationItem} from 'material-ui/BottomNavigation';
import Paper from 'material-ui/Paper';
import IconLocationOn from 'material-ui/svg-icons/communication/location-on';

const recentsIcon = <FontIcon className="material-icons">restore</FontIcon>;
const favoritesIcon = <FontIcon className="material-icons">favorite</FontIcon>;
const nearbyIcon = <IconLocationOn />;


export interface BottomClickActionsProps {  };
export interface BottomClickActionsState { selected: number };

export default class BottomClickActions extends React.Component<BottomClickActionsProps, BottomClickActionsState> {

    /**
     * Class constructor.
     */
    constructor() {
        super({});

    }

    componentWillMount() {
        this.setState({selected: 1});
    }

    select(bottomSelect: number)
    {
        this.setState({selected: bottomSelect});
    }

    /**
     * Render the component.
     */
    render() {


        return (

            <Paper zDepth={1}>
                    <BottomNavigation selectedIndex={this.state.selected}>
                    <BottomNavigationItem
                        label="Main"
                        icon={recentsIcon}
                        onClick={() => this.select(0)}
                    />
                    <BottomNavigationItem
                        label="Beispiel"
                        icon={favoritesIcon}
                        onClick={() => this.select(1)}
                    />
                    <BottomNavigationItem
                        label="Impressum/Rechtliche Informationen"
                        icon={nearbyIcon}
                        onClick={() => this.select(2)}
                    />
                    </BottomNavigation>
            </Paper>

        );
    }

}