import * as React from 'react'; 
import IconMenu from 'material-ui/IconMenu';
import IconButton from 'material-ui/IconButton';
import FontIcon from 'material-ui/FontIcon';
import NavigationExpandMoreIcon from 'material-ui/svg-icons/navigation/expand-more';
import MenuItem from 'material-ui/MenuItem';
import DropDownMenu from 'material-ui/DropDownMenu';
import RaisedButton from 'material-ui/RaisedButton';
import {Toolbar, ToolbarGroup, ToolbarSeparator, ToolbarTitle} from 'material-ui/Toolbar';


export interface MenuToolBarProps { currentSelection: any; locationLinks: any, onChange: any };
export interface MenuToolBarState { value: any; pageIdx: any; self: any; mode:any; };

export default class MenuToolBar extends React.Component<MenuToolBarProps,MenuToolBarState> {

  constructor(props) {
    super(props);


  }

  handleChange(event, index, value){


    this.setState({value: value});

    if (this.props.onChange)
    {
      this.props.onChange(value);
    }

  }


  componentWillMount()
  {
    this.setState({
      value: 0
    });
  }

  render() {

    var allPages = []

    for (var i = 0; i < this.props.locationLinks.length; i++)
    {
      allPages.push(<MenuItem key={i} value={this.props.locationLinks[i].index} primaryText={this.props.locationLinks[i].text}/>);
    }

    return (
      <Toolbar>
        <ToolbarGroup firstChild={true}>
          <DropDownMenu value={this.props.currentSelection} onChange={this.handleChange.bind(this)}>
            {allPages}
          </DropDownMenu>
        </ToolbarGroup>
      </Toolbar>
    );
  }
}