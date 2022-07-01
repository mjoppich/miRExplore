import * as React from 'react'; 
import AutoComplete from 'material-ui/AutoComplete';

import axios from 'axios';

import config from '../config';

/**
 *
 * `AutoComplete` search text can be implemented as a controlled value,
 * where `searchText` is handled by state in the parent component.
 * This value is reset with the `onNewRequest` callback.
 */
export interface AutoCompleteProps {onElementSelected: any};
export interface AutoCompleteState { searchText: string, terms: Array<string>};


export default class ACInput extends React.Component<AutoCompleteProps,AutoCompleteState> {
  
    constructor(props)
    {
        super(props);
    }

    componentWillMount()
    {
        this.setState({searchText: "", terms: []});
    }

  handleUpdateInput(searchText)
  {

    console.log("Search Text inserted: " + searchText)
    var self = this;

    self.setState({searchText: searchText})

    axios.post(config.getRestAddress() + "/autocomplete", {search: searchText}, config.axiosConfig)
          .then(function (response) {


            console.log(response.data)

            self.setState({terms: [].concat(response.data.genes).concat(response.data.mirna)})

          })
          .catch(function (error) {
            console.log(error)
            self.setState({terms: []})
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
          dataSource={this.state.terms}
          filter={(searchText, key) => (key.indexOf(searchText) !== -1)}
          openOnFocus={true}
        />
      </div>
    );
  }
}

