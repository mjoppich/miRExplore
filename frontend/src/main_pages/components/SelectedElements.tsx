import * as React from 'react'; 
import Avatar from 'material-ui/Avatar';
import Chip from 'material-ui/Chip';
import FontIcon from 'material-ui/FontIcon';
import SvgIconFace from 'material-ui/svg-icons/action/face';
import {blue300, indigo900} from 'material-ui/styles/colors';

const FaceIcon = <SvgIconFace/>



function handleRequestDelete() {
  alert('You clicked the delete button.');
}

function handleClick() {
  alert('You clicked the Chip.');
}

export interface SelectedElementsProps { elements: any };
export interface SelectedElementsState { };

const styles = {
  chip: {
    margin: 4,
  }};

/**
 * Examples of Chips, using an image [Avatar](/#/components/font-icon), [Font Icon](/#/components/font-icon) Avatar,
 * [SVG Icon](/#/components/svg-icon) Avatar, "Letter" (string) Avatar, and with custom colors.
 *
 * Chips with the `onRequestDelete` property defined will display a delete icon.
 */
export default class SelectedElements extends React.Component<SelectedElementsProps, SelectedElementsState> {

  render() {

    var allChips: Array<any> = [];

    for (var i = 0; i < this.props.elements.length; ++i)
    {
      let elem = this.props.elements[i];

      allChips.push(        <Chip
      key={i}
        backgroundColor={blue300}
        onRequestDelete={handleRequestDelete}
        onClick={handleClick}
        style={styles.chip}
      >
        <Avatar icon={FaceIcon} />
        {elem}
      </Chip>
      );
    }

    return (
      <div style={{display: 'flex', flexWrap: 'wrap'}}>

        <Chip
          onRequestDelete={handleRequestDelete}
          onClick={handleClick}
          style={styles.chip}
        >
          <Avatar icon={FaceIcon} />
          Deletable Avatar Chip
        </Chip>
      {allChips}
      </div>
    );
  }
}