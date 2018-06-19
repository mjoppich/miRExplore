import * as React from "react"; 

export interface ScrollAnchorProps { id: string }

export class ScrollAnchor extends React.Component<ScrollAnchorProps, {}> {

    constructor(props: ScrollAnchorProps) {
        super(props);

        console.log("Scroll ANchor");
        console.log(this.props);
    }

    componentDidMount()
    {
        console.log("Component Mounted with id " + this.props.id);
    }

    render() {
        return (
        <div>
        <div id={this.props.id} style={{position: "relative", top: "-80px"}}></div>
        <div>
            {this.props.children}
        </div>
        </div>);
    }
}