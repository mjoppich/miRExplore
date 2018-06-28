import * as React from 'react'; 

import ChipInput from 'material-ui-chip-input'
import AutoComplete from 'material-ui/AutoComplete'

import axios from 'axios';
import config from '../config';



export interface OrganismChipACProps { onValueChange: any};
export interface OrganismChipACState { allowedOrgs: Array<string>, values: Array<string>};


export default class OrganismChipAC extends React.Component<OrganismChipACProps, OrganismChipACState>{

    neo4jd3: any = null;

    constructor(props)
    {
        super(props);
    }

    componentWillMount()
    {
        var self=this;

        this.setState({allowedOrgs: []});

        axios.get(config.getRestAddress() + "/organisms", config.axiosConfig)
        .then(function (response) {
          console.log(response.data)

          self.setState({allowedOrgs: [].concat(response.data)})

        })
        .catch(function (error) {
          console.log(error)
          self.setState({allowedOrgs: []})
        });
    }

    componentDidMount(){

    }

    acceptOrganism(value)
    {
        if (this.state.allowedOrgs.indexOf(value) >= 0)
        {
            return true;
        }

        return false;
    }

    addOrganism(value)
    {
        this.state.values.push(value);
        this.props.onValueChange(this.state.values);
        this.setState({values: this.state.values})

    }
    deleteOrganism( elementText, i )
    {
        var idx = this.state.values.indexOf(elementText);

        if (idx >= 0)
        {
            this.state.values.splice(idx, 1);
        }
        this.props.onValueChange(this.state.values);
        this.setState({values: this.state.values})
    }

    handleOrganismAC(searchText)
    {

    }
  
    render(){
      return (

                <ChipInput
                            value={this.state.values}
                            onBeforeRequestAdd={(chip) => this.acceptOrganism(chip)}
                            onRequestAdd={(chip) => this.addOrganism(chip)}
                            onRequestDelete={(chip, index) => this.deleteOrganism(chip, index)}
                            openOnFocus={true}
                            filter={AutoComplete.fuzzyFilter}
                            dataSource={this.state.allowedOrgs}
                            onUpdateInput={(searchText, dataSource, params) => this.handleOrganismAC(searchText)}
                            fullWidth
                            fullWidthInput
                            floatingLabelText='Organism'
                            hintText='Organism Name'
                        />

        )
    }
}