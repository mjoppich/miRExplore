import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';

import getMuiTheme from 'material-ui/styles/getMuiTheme';
import MuiThemeProvider from 'material-ui/styles/MuiThemeProvider';

import MainPageBase from '../main_pages/container/MainPageBase';
import {WelcomePage as MainWelcomePage} from '../main_pages/container/MainPage';

//@ts-ignore
import ScrollableAnchor from 'react-scrollable-anchor';

import FontIcon from 'material-ui/FontIcon';
import {BottomNavigation, BottomNavigationItem} from 'material-ui/BottomNavigation';
import Paper from 'material-ui/Paper';

import ActionExplore from 'material-ui/svg-icons/action/explore';
import ActionHome from 'material-ui/svg-icons/action/home';
import ActionTimeline from 'material-ui/svg-icons/action/timeline';
import MenuToolBar from '../main_pages/components/MenuToolBar';


const ExploreIcon = <ActionExplore/>
const HomeIcon = <ActionHome/>
const AnalyzeIcon = <ActionTimeline/>


export interface MainPageProps { selected: number};
export interface MainPageState { selected: number, subPageElems: any, locLinks: any, width: number, height: number };

export default class MainPage extends React.Component<MainPageProps, MainPageState> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);


    }

    componentWillMount() {

        var mainLocLinks = [];
        mainLocLinks.push( {index: 0, text: 'Home'} );

        const mainPage = (
            <div>        
                <div style={{margin: "15px auto 15px auto"}}><MainWelcomePage switchTab={this.select.bind(this)}/></div>
            </div>
        );


        var allLinks = [];
        mainLocLinks.forEach(element => { allLinks.push(element); });


        this.setState({
            subPageElems: [
                mainPage
            ],
            locLinks: allLinks,
            selected: this.props.selected
        })
    }

    select(bottomSelect: number)
    {

        if ((bottomSelect != null) && (bottomSelect != undefined))
        {
            this.setState({selected: bottomSelect});
        }
    }

    componentDidMount() {

        this.updateWindowDimensions();
        window.addEventListener('resize', this.updateWindowDimensions.bind(this));
    }

    componentWillUnmount() {
        window.removeEventListener('resize', this.updateWindowDimensions.bind(this));
    }

    updateWindowDimensions() {
        this.setState({ width: window.innerWidth, height: window.innerHeight });
    }

    /**
     * Render the component.
     */
    render() {
        var self=this;

        var expContainerWidth: number = 1000;
        var containerWidth = expContainerWidth + 'px';

        if (this.state.width*0.8 < expContainerWidth)
        {
            containerWidth = "80%"
        }

        containerWidth = "80%";

        return (
        <MuiThemeProvider muiTheme={getMuiTheme()}>
            <div>
                <div style={{height: 100, left:0, top: 0, width: 100+'%', zIndex: 1000, position: 'fixed'}}>
                    <MenuToolBar locationLinks={this.state.locLinks} currentSelection={this.state.selected} onChange={self.select.bind(self)}/>
                </div>

                <div style={{height: 100, top: 56+'px', width: containerWidth, zIndex: 900, position: 'relative', margin: '0 auto'}}>
                    {this.state.subPageElems[this.state.selected]}

                    <Card>
                            <BottomNavigation selectedIndex={this.state.selected}>
                            <BottomNavigationItem
                                label="Home"
                                icon={HomeIcon}
                                onClick={() => this.select(0)}
                            />

                            </BottomNavigation>
                    </Card>
                </div>

            </div>
        </MuiThemeProvider>);
      }

}