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
import { ExploreMainPage } from "main_pages/container/ExploreMainPage";
import {AnalyseMainPage} from "main_pages/container/AnalyseMainPage";

const ExploreIcon = <ActionExplore/>
const HomeIcon = <ActionHome/>
const AnalyzeIcon = <ActionTimeline/>


export interface MainPageProps { selected:number };
export interface MainPageState { selected: number, subPageElems: any; locLinks };

export default class MainPage extends React.Component<MainPageProps, MainPageState> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

        console.log(this.props);

    }

    componentWillMount() {

        var mainLocLinks = [];
        mainLocLinks.push( {href: '#welcome', text: 'Welcome !'} );

        var exampleLocLinks = [];
        exampleLocLinks.push( {href: '#explore', text: 'Welcome !'} );

        var rightsLocLinks = [];
        rightsLocLinks.push( {href: '#welcome', text: 'Welcome !'} );


        const mainPage = (
            <div>        
                <ScrollableAnchor id={'welcome'}><div style={{margin: "15px auto 15px auto"}}><MainWelcomePage/></div></ScrollableAnchor>
            </div>
        );

        const examplePage = (
            <div>
                <ScrollableAnchor id={'welcome'}><div style={{margin: "15px auto 15px auto"}}><ExploreMainPage/></div></ScrollableAnchor>
            </div>
        );

        const rightsPage = (
            <div>        
                <ScrollableAnchor id={'welcome'}><div style={{margin: "15px auto 15px auto"}}><AnalyseMainPage/></div></ScrollableAnchor>
            </div>
        );

        this.setState({
            subPageElems: [
                mainPage, examplePage, rightsPage
            ],
            locLinks: [
                mainLocLinks, exampleLocLinks, rightsLocLinks
            ],
            selected: this.props.selected
        })
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
        <MuiThemeProvider muiTheme={getMuiTheme()}>

            <MainPageBase locLinks={this.state.locLinks[this.state.selected]}>
                {this.state.subPageElems[this.state.selected]}
                <Card>
                        <BottomNavigation selectedIndex={this.state.selected}>
                        <BottomNavigationItem
                            label="Home"
                            icon={HomeIcon}
                            onClick={() => this.select(0)}
                        />
                        <BottomNavigationItem
                            label="Explore"
                            icon={ExploreIcon}
                            onClick={() => this.select(1)}
                        />
                        <BottomNavigationItem
                            label="Analyze"
                            icon={AnalyzeIcon}
                            onClick={() => this.select(2)}
                        />
                        </BottomNavigation>
                </Card>
            </MainPageBase>
        </MuiThemeProvider>);
      }

}