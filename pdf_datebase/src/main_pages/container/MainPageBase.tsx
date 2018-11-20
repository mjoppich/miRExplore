import * as React from "react"; 
import { Link } from 'react-router-dom';

import MenuToolBar from '../components/MenuToolBar';

export interface PageBaseProps { children: any; locLinks: any; };
export interface PageBaseState { width: any; height: any; hashLocation: any; self: any; mode:any; children: any; locLinks: any; };

export default class MainPageBase extends React.Component<PageBaseProps,PageBaseState> {

    constructor(props) {
    super(props);

    console.log(props);

    this.updateWindowDimensions = this.updateWindowDimensions.bind(this);
    }

    componentWillMount()
    {
        this.setState({ width: '0', height: '0', hashLocation: '', children: this.props.children, locLinks: this.props.locLinks });

    }

    componentDidMount() {

        this.updateWindowDimensions();
        window.addEventListener('resize', this.updateWindowDimensions);
    }

    componentWillUnmount() {
    window.removeEventListener('resize', this.updateWindowDimensions);
    }

    updateWindowDimensions() {

        //console.log("inner width: " + window.innerWidth);

    this.setState({ width: window.innerWidth, height: window.innerHeight });
    }

  componentWillReceiveProps(nextProps){
      //console.log("pagebase receive Props");
        console.log(nextProps);
        console.log(window.location);

        this.setState(nextProps);
        this.setState({hashLocation: window.location.hash});

  }

    render() {

        var expContainerWidth: number = 1000;
        var containerWidth = expContainerWidth + 'px';

        if (this.state.width*0.8 < expContainerWidth)
        {
            containerWidth = "80%"
        }

        containerWidth = "80%";

        return (
<div style={{display: 'flex', justifyContent: 'center'}}>

                <div style={{width: 0.8*this.state.width}}>
                    <div style={{height: 100, left:0, top: 0, width: 100+'%', zIndex: 1000, position: 'fixed'}}>
                        <MenuToolBar locationLinks={this.state.locLinks} currentSelection={this.state.hashLocation} onChange={null}/>
                    </div>

                    <div style={{height: 100, top: 56+'px', width: containerWidth, zIndex: 900, position: 'relative', margin: '0 auto'}}>
                        {this.state.children}
                    </div>

                </div>
            </div>
        );
    }
};
