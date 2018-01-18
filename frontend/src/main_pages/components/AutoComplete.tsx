import * as React from 'react'; 
import AutoComplete from 'material-ui/AutoComplete';

const colors = [
  'Red',
  'Orange',
  'Yellow',
  'Green',
  'Blue',
  'Purple',
  'Black',
  'White',
];

/**
 * `AutoComplete` search text can be implemented as a controlled value,
 * where `searchText` is handled by state in the parent component.
 * This value is reset with the `onNewRequest` callback.
 */
export interface AutoCompleteProps {onElementSelected: any};
export interface AutoCompleteState { searchText: string};


export default class ACInput extends React.Component<AutoCompleteProps,AutoCompleteState> {
  
    constructor(props)
    {
        super(props);
    }

    componentWillMount()
    {
        this.setState({searchText: ""});
    }

  handleUpdateInput(searchText)
  {
    this.setState({
      searchText: searchText,
    });
  };

  handleNewRequest(chosenRequest: string, index: number) {
    this.setState({
      searchText: '',
    });

    // fire onUpdateEvent
    this.props.onElementSelected(chosenRequest);
  };

  render() {
    return (
      <div>
        <AutoComplete
          hintText="Type 'r', case insensitive"
          searchText={this.state.searchText}
          onUpdateInput={this.handleUpdateInput.bind(this)}
          onNewRequest={this.handleNewRequest.bind(this)}
          dataSource={colors}
          filter={(searchText, key) => (key.indexOf(searchText) !== -1)}
          openOnFocus={true}
        />
      </div>
    );
  }
}

