import * as React from "react"; 
import { Card, CardTitle, CardText } from 'material-ui/Card';
import Dropzone from 'react-dropzone';
import ReactDOMServer from 'react-dom/server';
import FlatButton from 'material-ui/FlatButton';
import axios from 'axios';
import config from '../config';
import ReactTable from 'react-table';
import matchSorter from 'match-sorter';
import Button from '@material-ui/core/Button';
import LinearProgress from '@material-ui/core/LinearProgress';
import D3SVGParallelLinesGraph from '../components/D3SVGForceParallelLines';

import OboChipAC from '../components/OBOChipAC';

import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel'

export interface WelcomePageState { files: any, results: any, authorResults: any, queryState:any, interactionData: any, updateInteractions: boolean, selectedAuthors: Array<any>, }
export class WelcomePage extends React.Component<{ switchTab?: any },WelcomePageState> {

    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);

    }

    /**
     * This method will be executed after initial rendering.
     */
    componentDidMount() {

    }

    resid2file = {};



    readonly state = {
        'files': [],
        'results': [],
        'queryState': 0,
        interactionData: [],
        updateInteractions: false,
        selectedAuthors: [],
        authorResults: {}
    };
    
      onDrop(files) {
        this.setState({
          files: files
        });
      }

      sendFilesToServer()
      {

        if (this.state.queryState > 0)
        {
            return;
        }

        var formData = new FormData();
        this.resid2file = {};

        for (var i = 0; i < this.state.files.length; ++i)
        {
            var elem = this.state.files[i];

            var elemid = "pdf"+i

            this.resid2file[elemid] = this.state.files[i].name
            console.log(elemid);
            formData.append(elemid, elem);

        }

        var self=this;
        this.setState({queryState: 1})

        axios.post(config.getRestAddress() + "/textmine",formData, config.axiosConfig)
        .then(function (response) {

            console.log("Received TM Results")
            console.log(response.data)

            self.setState({results: response.data, queryState: 0, updateInteractions: true});
        })
        .catch(function (error) {
          console.log(error);
          self.setState({results: [], queryState: -1, updateInteractions: false})
        });


    }

    makeTMResult()
    {



        if (this.state.queryState == 1)
        {
            var queryOngoing = [];
            queryOngoing.push(<LinearProgress key={queryOngoing.length}/>);
            queryOngoing.push(<p key={queryOngoing.length}>We are currently processing your PDF files.</p>);

            return (<Card style={{marginBottom: "20px"}}>
                    <CardText>
                        {queryOngoing}
                    </CardText>
            </Card>);

        } else if (this.state.queryState < 0)
        {
            var queryOngoing = [];
            queryOngoing.push(<p key={queryOngoing.length}>An error occurred processing your files. Please try again.</p>);

            return (<Card style={{marginBottom: "20px"}}>
                    <CardText>
                        {queryOngoing}
                    </CardText>
            </Card>);
        }

        if (this.state.results.length == 0)
        {
            return null;
        }

        var allPDFTables = [];

        var allResPDFs = Object.keys(this.state.results);

        for (var i = 0; i < allResPDFs.length; ++i)
        {

            var pdfData = this.state.results[allResPDFs[i]];

            allPDFTables.push(<h2 key={allPDFTables.length}>{"Gene-mRNA interactions: " + this.resid2file[allResPDFs[i]]}</h2>)

                allPDFTables.push(<ReactTable
                    key={allPDFTables.length}
                    data={pdfData}
                    filterable
                    defaultFilterMethod={(filter, row) =>
                      String(row[filter.id]) === filter.value}
                    columns={[
                      {
                        Header: "Found Interactions",
                        columns: [
                            {
                                Header: "Gene ID",
                                id: "gene",
                                accessor: d => d.gene,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['gene'] ];
                                    
                                    var retval = matchSorter(elems, filter.value);
        
                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    //var sects = row.value.split(".").slice(-2)

                                    return <span>{row.value}</span>
                                }
                            },
                            {
                                Header: "miRNA ID",
                                id: "mirna",
                                accessor: d => d.mirna,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [rows['mirna']];
                                    
                                    var retval = matchSorter(elems, filter.value);
        
                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    return <span>{row.value}</span>
                                }
                            },
                        ]},
                        { Header: "Evidences",
                        columns: [
                            {
                                Header: "Pubmed Evidences",
                                id: "pub_evs",
                                accessor: d => d.pubevs,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [];

                                    if (rows["pubevs"])
                                    {
                                        rows['pubevs'].forEach(element => {
                                            elems.push(element[0])

                                            element[3].forEach(elem => {

                                                var author = elem[0] + ", " + elem[2];
    
                                                if (elem[1] && elem[1].length > 0)
                                                {
                                                    author += " " + elem[1] 
                                                }
    
                                                elems.push(author)
                                            })
                                        });

                                    }

                                    var retval = matchSorter(elems, filter.value);

                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    var allDisInfo = [];
                                
                                    var relInfo = row.value;

                                    relInfo.sort(function(a,b) {return a[0] > b[0];})
    
                                    //console.log(relInfo);    
                                    for (var i = 0; i < relInfo.length; ++i)
                                    {   

                                        var linkInfo = relInfo[i];

                                        var authors = [];

                                        linkInfo[3].forEach(elem => {

                                            var author = elem[0] + ", " + elem[2];

                                            if (elem[1] && elem[1].length > 0)
                                            {
                                                author += " " + elem[1] 
                                            }

                                            authors.push(author)
                                        })

                                        allDisInfo.push(
                                                <div>
                                                    <span key={i} style={{display: "block", whiteSpace: "normal"}}>                                        
                                                        <a href={"https://www.ncbi.nlm.nih.gov/pubmed/"+linkInfo[0]}>{linkInfo[2] + " (" + linkInfo[1] + ")"}</a>
                                                    </span><br/>
                                                    <span key={i} style={{display: "block", whiteSpace: "normal"}}>                                        
                                                        {authors.join(" and ")}
                                                    </span><br/>
                                                </div>
                                            );

                                        
                                    }
    
                                    return <div>{allDisInfo}</div>;
                                }
                            },
                            {
                                Header: "Gene Text Context",
                                id: "gene_context",
                                accessor: d => d.gene_context,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['gene_context'] ];
                                    
                                    var retval = matchSorter(elems, filter.value);

                                    return retval.length > 0;
                                },
                                Cell: (row) => {
                                    return <span style={{display: "block", whiteSpace: "normal"}}>{row.value}</span>
                                }
                            },
                            {
                                Header: "miRNA Text Context",
                                id: "mirna_context",
                                accessor: d => d.mirna_context,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['mirna_context'] ];
                                    
                                    var retval = matchSorter(elems, filter.value);

                                    return retval.length > 0;
                                },
                                Cell: (row) => {
                                    return <span style={{display: "block", whiteSpace: "normal"}}>{row.value}</span>
                                }
                            },
                            
                        ]
                      }
                    ]}
                    defaultPageSize={10}
                    className="-striped -highlight"
                    
                  />);



            console.log("Spongebob")
            console.log(this.state);

        }



        
        return (<div>
                    <Card style={{marginBottom: "20px"}}>
                        <CardText>
                            {allPDFTables}
                        </CardText>
                    </Card>
                </div>
                
        );

    }


    queryAuthors()
    {

      if (this.state.queryState > 0)
      {
          return;
      }

      var self=this;
      this.setState({queryState: 1})

      axios.post(config.getRestAddress() + "/authorsearch", {authors: self.state.selectedAuthors}, config.axiosConfig)
      .then(function (response) {

          console.log("Received TM Results")
          console.log(response.data)

          self.setState({authorResults: response.data, queryState: 0, updateInteractions: true});
      })
      .catch(function (error) {
        console.log(error);
        self.setState({results: [], queryState: -1, updateInteractions: false})
      });


  }

    makeAuthorResult()
    {

        if (Object.keys(this.state.authorResults).length == 0)
        {
            return <div></div>;
        }

        var author2rels = this.state.authorResults["author2rels"];
        console.log(author2rels);

        if (author2rels == undefined)
        {
            return <div></div>;
        }

        var relKeys = Object.keys(author2rels);

        var allResElems = [];

        for (var i = 0; i < relKeys.length; ++i)
        {
            allResElems.push(<h2 key={allResElems.length}>{relKeys[i]}</h2>)

            var relEvs = author2rels[relKeys[i]];

            allResElems.push(<ReactTable
                key={allResElems.length}
                data={relEvs}
                filterable
                defaultFilterMethod={(filter, row) =>
                  String(row[filter.id]) === filter.value}
                columns={[
                  {
                    Header: "Found Interactions",
                    columns: [
                        {
                            Header: "Gene ID",
                            id: "gene",
                            accessor: d => d.gene,
                            filterMethod: (filter, rows) => 
                            {
                                var elems = [ rows['gene'] ];
                                
                                var retval = matchSorter(elems, filter.value);
    
                                return retval.length > 0;
                            },
                            Cell: (row) => {

                                //var sects = row.value.split(".").slice(-2)

                                return <span>{row.value}</span>
                            }
                        },
                        {
                            Header: "miRNA ID",
                            id: "mirna",
                            accessor: d => d.mirna,
                            filterMethod: (filter, rows) => 
                            {
                                var elems = [rows['mirna']];
                                
                                var retval = matchSorter(elems, filter.value);
    
                                return retval.length > 0;
                            },
                            Cell: (row) => {

                                return <span>{row.value}</span>
                            }
                        },
                    ]},
                    { Header: "Evidences",
                    columns: [
                        {
                            Header: "Pubmed Evidences",
                            id: "pub_evs",
                            accessor: d => d.pubevs,
                            filterMethod: (filter, rows) => 
                            {
                                var elems = [];

                                if (rows["pubevs"])
                                {
                                    rows['pubevs'].forEach(element => {
                                        elems.push(element[0])

                                        element[3].forEach(elem => {

                                            var author = elem[0] + ", " + elem[2];

                                            if (elem[1] && elem[1].length > 0)
                                            {
                                                author += " " + elem[1] 
                                            }

                                            elems.push(author)
                                        })
                                    });

                                }

                                var retval = matchSorter(elems, filter.value);

                                return retval.length > 0;
                            },
                            Cell: (row) => {

                                var allDisInfo = [];
                            
                                var relInfo = row.value;

                                relInfo.sort(function(a,b) {return a[0] > b[0];})

                                //console.log(relInfo);    
                                for (var i = 0; i < relInfo.length; ++i)
                                {   

                                    var linkInfo = relInfo[i];

                                    var authors = [];

                                    linkInfo[3].forEach(elem => {

                                        var author = elem[0] + ", " + elem[2];

                                        if (elem[1] && elem[1].length > 0)
                                        {
                                            author += " " + elem[1] 
                                        }

                                        authors.push(author)
                                    })

                                    allDisInfo.push(
                                            <div>
                                                <span key={i} style={{display: "block", whiteSpace: "normal"}}>                                        
                                                    <a href={"https://www.ncbi.nlm.nih.gov/pubmed/"+linkInfo[0]}>{linkInfo[2] + " (" + linkInfo[1] + ")"}</a>
                                                </span><br/>
                                                <span key={i} style={{display: "block", whiteSpace: "normal"}}>                                        
                                                    {authors.join(" and ")}
                                                </span><br/>
                                            </div>
                                        );

                                    
                                }

                                return <div>{allDisInfo}</div>;
                            }
                        },
                        {
                            Header: "Gene Text Context",
                            id: "gene_context",
                            accessor: d => d.gene_context,
                            filterMethod: (filter, rows) => 
                            {
                                var elems = [ rows['gene_context'] ];
                                
                                var retval = matchSorter(elems, filter.value);

                                return retval.length > 0;
                            },
                            Cell: (row) => {
                                return <span style={{display: "block", whiteSpace: "normal"}}>{row.value}</span>
                            }
                        },
                        {
                            Header: "miRNA Text Context",
                            id: "mirna_context",
                            accessor: d => d.mirna_context,
                            filterMethod: (filter, rows) => 
                            {
                                var elems = [ rows['mirna_context'] ];
                                
                                var retval = matchSorter(elems, filter.value);

                                return retval.length > 0;
                            },
                            Cell: (row) => {
                                return <span style={{display: "block", whiteSpace: "normal"}}>{row.value}</span>
                            }
                        },
                        
                    ]
                  }
                ]}
                defaultPageSize={10}
                className="-striped -highlight"
                
              />);

        }


        return <div>{allResElems}</div>;
    }

    /**
     * Render the component.
     */
    render() {

        var tmresult = this.makeTMResult();
        var authorResult = this.makeAuthorResult();

        const dropzoneStyle = {
            width  : "100%",
            height : "100px",
            border : "1px dashed black"
        };

        return (

            <div>
            <Card style={{marginBottom: "20px"}}>
                <CardTitle
                    title="miRExplore DateBase - Welcome!"
                    subtitle="Want to know who first discovered a miRNA-Gene interaction?"
                />
                <CardText >
                    <p>Drop your PDFs and scan for miRNA-Target Interactions - and find the original authors.</p>
                </CardText>
            </Card>
            <Card style={{marginBottom: "20px", backgroundColor: "#ff786e"}}>
                <CardTitle
                    title="PDF Drag'n'Drop Zone"
                    subtitle=""
                />
                <CardText>
                    <p>
                        Drop your PDFs and scan for miRNA-Target Interactions.
                    </p>
                    <p>
                        You may drop multiple PDFs at one time!
                    </p>
                </CardText>
            </Card>
            <Card style={{marginBottom: "20px"}}>
                <CardTitle>DropZone for your PDFs</CardTitle>
                <CardText >
                <div>
                    <section>
                        <div className="dropzone">
                        <Dropzone onDrop={this.onDrop.bind(this)} accept="application/pdf" style={dropzoneStyle}>
                            <p>Try dropping some files here, or click to select files to upload.</p>
                        </Dropzone>
                        </div>
                    </section>

                    <p>Considered Files:</p>
                    <section>
                        <ul>
                            {
                            this.state.files.map(f => <li key={f.name}>{f.name} - {f.size} bytes</li>)
                            }
                        </ul>
                    </section>

                    <FormGroup>
                    <Button onClick={() => this.sendFilesToServer()}>
                    Submit To Server
                    </Button>
                    </FormGroup>


                </div>


                </CardText>
            </Card>
            {tmresult}

           <Card style={{marginBottom: "20px"}}>
                <CardTitle>Search For Author</CardTitle>
                <CardText >            

                        <OboChipAC
                            url="authorcomplete"
                            floatText="Find Interactions for Authors"
                            hintText="Author First/Last Name"
                            onValueChange={(newvalues) => this.setState({selectedAuthors: newvalues})
                        }/>

                        <Button onClick={() => this.queryAuthors()}>
                            Query Authors
                        </Button>

                </CardText>
            </Card>
            {
                authorResult
            }

            <Card style={{marginBottom: "20px"}}>
                <CardTitle
                    title="Data Privacy"
                />
                <CardText >
                    <p>You can find the data privacy guide-lines <a target="_blank" href="https://www.bio.ifi.lmu.de/funktionen/datenschutz/index.html">here</a>.</p>
                    <p>Your IP address might be stored by the web server. All log files are deleted after 7 days.</p>
                </CardText>
            </Card>
            </div>

        );
    }

}