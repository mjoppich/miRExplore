
import axios from 'axios';
import config from '../config';

import Autocomplete from "@mui/material/Autocomplete";
import TextField from "@mui/material/TextField";
import Chip from "@mui/material/Chip";

import React, { useEffect, useState } from "react";

interface IEntityElement {
    name: string;
    group: string;
 }

export interface OboChipACProps { detailsurl:string, url: string, floatText:string, hintText: string, onValueChange: any};
export interface OboChipACState { values: Array<string>, options: Array<string>, elements: Array<IEntityElement>};


export default class OboChipAC extends React.Component<OboChipACProps, OboChipACState>{

    neo4jd3: any = null;

    constructor(props)
    {
        super(props);
        this.state = {values: [], options: [], elements:[]};

    }

    componentDidMount(){

    }


    addElement(values)
    {

        console.log("New List of Values")
        console.log(values)
        var self = this;

        axios.post(config.getRestAddress() + "/" + this.props.detailsurl, {search: values}, config.axiosConfig)
        .then(function (response) {
            console.log(response.data)

            var allElements = response.data;
            var elemnames = [];
            allElements.forEach(element => {
                elemnames.push(element.name)
            });

            self.setState({elements: allElements, values: elemnames})

            self.props.onValueChange(allElements);

        })
        .catch(function (error) {

            console.log(config.getRestAddress() + "/" + self.props.detailsurl)
            console.log(values)
            console.log(error)
        });

    }


    handleAutoComplete(searchText)
    {
        var self = this;

        console.log(searchText)

        axios.post(config.getRestAddress() + "/" + this.props.url, {search: searchText}, config.axiosConfig)
        .then(function (response) {
            console.log(response.data)

            var allElements = response.data;
            var elemnames = [];
            allElements.forEach(element => {
                elemnames.push(element.name)
            });

            self.setState({options: elemnames})

        })
        .catch(function (error) {

            console.log(config.getRestAddress() + "/" + self.props.url)
            console.log(searchText)
            console.log(error)
            self.setState({options: []})
        });
    }
  

    render(){
        var self=this;

        return (
            <Autocomplete
              multiple
              id="tags-filled"  
              options={self.state.options}
              defaultValue={[]}
              value={self.state.values}
              freeSolo
              onChange={(e:any, value:any) => self.addElement(value)}
              onInputChange={(e:any, value:any, r:any) => self.handleAutoComplete(value)}
              renderTags={(
                value: any[],
                getTagProps: (arg0: { index: any }) => JSX.IntrinsicAttributes
              ) =>
                value.map((option: any, index: any) => {
                  return (
                    <Chip
                      key={index}
                      variant="outlined"
                      label={option}
                      {...getTagProps({ index })}
                    />
                  );
                })
              }
              renderInput={(params: any) => (
                <TextField
                  {...params}
                  label={self.props.floatText}
                  placeholder="Add a receiver by pressing enter after its dotName or address"
                />
              )}
              style={{padding: "10px" }}
            />)
    }
}