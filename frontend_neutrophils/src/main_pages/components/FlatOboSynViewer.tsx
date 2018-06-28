import * as React from "react"; 
import axios from 'axios';
import config from '../config';
import ReactTable from 'react-table';
import matchSorter from 'match-sorter';


export interface FlatOboSynViewerProps { oboname: string }
export interface FlatOboSynViewerState { oboinfo: any }



export default class FlatOboSynViewer extends React.Component<FlatOboSynViewerProps, FlatOboSynViewerState> {

    constructor(props) {
        super(props);

        console.log(this.props);

    }


    readonly state = {
        'oboinfo': null
    };

    componentDidMount()
    {
        console.log("NI Sankey did mount");
        this.reloadGraphData();
    }

    componentDidUpdate(prevProps, prevState, snapshot) {

        if (prevProps != this.props)
        {
            this.reloadGraphData();
        }
    }


    reloadGraphData()
    {


        if (this.props.oboname == null)
        {
            return;
        }

        var sendData = {};
        var self=this;


        sendData['oboname'] = this.props.oboname;

        axios.post(config.getRestAddress() + "/oboinfo",sendData, config.axiosConfig)
        .then(function (response) {

            console.log("Oboinfo data fetch")
            console.log(response.data)

            self.setState({oboinfo: response.data});
            
        })
        .catch(function (error) {
          console.log(error);
          self.setState({oboinfo: null});
        });
    }

    render() {

        var content = <span>No information found for obo {this.props.oboname}</span>;

        if (this.state.oboinfo != null)
        {

            content =  <ReactTable
                    filterable
                    defaultFilterMethod={(filter, row) => String(row[filter.id]) === filter.value}
                    data={this.state.oboinfo}
                    columns={
                        [

                            {
                                Header: "Ontology",
                                columns: [
                                    {
                                        Header: "Term ID/Name",
                                        id:"termid",
                                        accessor: d => [d.termid, d.termname],
                                        Cell: (row) => {
                                            var docID = row.value[0];
                                            var docSent = row.value[1];

                                            var allHTMLElems = [];
                                            allHTMLElems.push(<a key={allHTMLElems.length} href={"https://www.ncbi.nlm.nih.gov/pubmed/"+docID} target="_blank">{docID}</a>);
                                            allHTMLElems.push(<span key={allHTMLElems.length}>{docSent}</span>);
                                            
                                            return <div style={{display: "grid"}}>
                                                {
                                                    allHTMLElems
                                                }
                                            </div>;
                                        },
                                        filterMethod: (filter, row) => {
                                            var filterID = filter.id;
                                            var rowData = row[filterID];

                                            var allTerms = [];
                                            allTerms.push(rowData[0]);
                                            allTerms.push(rowData[1]);
                                            
                                            var retval = matchSorter(allTerms, filter.value);
                                            return retval.length > 0;
                                        }
                                    },
                                    {
                                        Header: "Children",
                                        id:"children",
                                        accessor: d => d.children,
                                        Cell: (row) => {
                                            var children = row.value;    
                                            var allHTMLElems = [];

                                            for (var i =0; i < children.length; ++i)
                                            {
                                                allHTMLElems.push(<a key={allHTMLElems.length} href={"https://www.ncbi.nlm.nih.gov/pubmed/"+children[i].id} target="_blank">{children[i].name}</a>);
                                            }
                                            
                                            return <div style={{display: "grid"}}>
                                                {
                                                    allHTMLElems
                                                }
                                            </div>;
                                        },
                                        filterMethod: (filter, row) => {
                                            var filterID = filter.id;
                                            var rowData = row[filterID];

                                            var allTerms = [];
                                            for (var i =0; i < rowData.length; ++i)
                                            {
                                                allTerms.push(rowData[i].name)
                                            }
                                            
                                            var retval = matchSorter(allTerms, filter.value);
                                            return retval.length > 0;
                                        }
                                    },
                                    {
                                        Header: "Synonymes",
                                        id:"syns",
                                        accessor: d => d.synonymes,
                                        Cell: (row) => {
                                            var syns = row.value;    
                                            var allHTMLElems = [];

                                            for (var i =0; i < syns.length; ++i)
                                            {
                                                allHTMLElems.push(<span key={allHTMLElems.length}>{syns[i]}</span>);
                                            }
                                            
                                            return <div style={{display: "grid"}}>
                                                {
                                                    allHTMLElems
                                                }
                                            </div>;
                                        },
                                        filterMethod: (filter, row) => {
                                            var filterID = filter.id;
                                            var rowData = row[filterID];

                                            var allTerms = [];
                                            for (var i =0; i < rowData.length; ++i)
                                            {
                                                allTerms.push(rowData[i])
                                            }
                                            
                                            var retval = matchSorter(allTerms, filter.value);
                                            return retval.length > 0;
                                        }
                                    }


                                    
                                ]
                            }
                        ]
                    }
                    defaultPageSize={10}
                    showPaginationTop
                    showPaginationBottom
                    className="-striped -highlight"

                    />



        }


        return (<div>
            {content}
        </div>);
    }

}