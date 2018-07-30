import * as React from "react"; 
import { Router, Route, Switch } from 'react-router';


import createBrowserHistory from 'history/createBrowserHistory';


//@ts-ignore
import {configureAnchors} from 'react-scrollable-anchor';

/*
//@ts-ignore
import * as tep from 'react-tap-event-plugin';
*/

import MainPage from './MainPage';
import NotFoundPage from './NotFoundPage';


var injectTapEventPlugin = require("react-tap-event-plugin");
injectTapEventPlugin();
configureAnchors({offset: -60});

//var routerHistory = createBrowserHistory();

export class RoutingInformation
{
    exactPath:string;
    renderFunc: any;

    constructor(exactPath: string, renderFunc: any)
    {
        this.exactPath = exactPath;
        this.renderFunc = renderFunc;
    }

}
  
export class MainApp extends React.Component<{},{}> {
    render() {

        var allRouteInfos = [];
        allRouteInfos.push(new RoutingInformation("/", (props) => <MainPage selected={0} {...props} />))
        allRouteInfos.push(new RoutingInformation("/explore", (props) => <MainPage selected={1} {...props} />))
        allRouteInfos.push(new RoutingInformation("/analyse", (props) => <MainPage selected={2} {...props} />))

        var allRoutes = []

        for (var i = 0; i < allRouteInfos.length; ++i)
        {

            var routeInfo : RoutingInformation = allRouteInfos[i];

            allRoutes.push(<Route key={i} exact path={routeInfo.exactPath} render={routeInfo.renderFunc} />)
        }

        return (
        <Switch>
            {allRoutes}
            <Route component={NotFoundPage} />
        </Switch>
        );
    }
}