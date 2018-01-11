import * as React from 'react'; 
import IconMenu from 'material-ui/IconMenu';
import IconButton from 'material-ui/IconButton';
import FontIcon from 'material-ui/FontIcon';
import NavigationExpandMoreIcon from 'material-ui/svg-icons/navigation/expand-more';
import MenuItem from 'material-ui/MenuItem';
import DropDownMenu from 'material-ui/DropDownMenu';
import RaisedButton from 'material-ui/RaisedButton';
import {Toolbar, ToolbarGroup, ToolbarSeparator, ToolbarTitle} from 'material-ui/Toolbar';


export interface MenuToolBarProps {  currentLocation: any; locationLinks: any };
export interface MenuToolBarState { value: any; pages: any; pageIdx: any; self: any; mode:any; };

export default class MenuToolBar extends React.Component<MenuToolBarProps,MenuToolBarState> {

  constructor(props) {
    super(props);


  }

  handleChange(event, index, value){
       this.setState({value});

  }

  componentWillReceiveProps(nextProps){

      if ('locationLinks' in nextProps)
      {
        var href2idx = {}
        for (var i = 0; i < nextProps.locationLinks.length; i++)
        {
          href2idx[ nextProps.locationLinks[i].href ] = i;
        }

        this.setState({
          pages: nextProps.locationLinks,
          pageIdx: href2idx,
          value: 0
        });
    
        console.log("Updated locationLinks in MenuToolBar");

      } else {
        var newlocation = nextProps.currentLocation || ''

        if (newlocation === '')
        {
          return;
        }
  
        var idx = this.state.pageIdx[ newlocation ] || 0;
        this.setState({value: idx});
      }



      //console.log("toolbar " + idx + " " + newlocation);

  }

  componentWillMount()
  {
    var href2idx = {}
    for (var i = 0; i < this.props.locationLinks.length; i++)
    {
      href2idx[ this.props.locationLinks[i].href ] = i;
    }

    this.setState({
      value: 0,
      pages: this.props.locationLinks,
      pageIdx: href2idx
    });
  }

  render() {

    var allPages = []

    for (var i = 0; i < this.state.pages.length; i++)
    {
      allPages.push(<MenuItem key={i} value={i} primaryText={this.state.pages[i].text} href={this.state.pages[i].href}/>);
    }

    return (
      <Toolbar>
        <ToolbarGroup firstChild={true}>
          <DropDownMenu value={this.state.value} onChange={this.handleChange.bind(this)}>
            {allPages}
          </DropDownMenu>
        </ToolbarGroup>
      </Toolbar>
    );
  }
}