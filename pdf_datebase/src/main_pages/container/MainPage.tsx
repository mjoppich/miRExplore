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

import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel'

export interface WelcomePageState { files: any, results: any, queryState:any, interactionData: any, updateInteractions: boolean }
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
        updateInteractions: false
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

            var pdfConcepts = Object.keys(pdfData);

            for (var j = 0; j < pdfConcepts.length; ++j)
            {

                var contextData = pdfData[pdfConcepts[j]];

                if ((contextData.length ==  0) || (contextData.results.length == 0))
                {
                    continue;
                }

                allPDFTables.push(<h2 key={allPDFTables.length}>{this.resid2file[allResPDFs[i]]} - {pdfConcepts[j]}</h2>)

                allPDFTables.push(<ReactTable
                    key={allPDFTables.length}
                    data={contextData.results}
                    filterable
                    defaultFilterMethod={(filter, row) =>
                      String(row[filter.id]) === filter.value}
                    columns={[
                      {
                        Header: "Found Concepts",
                        columns: [
                            {
                                Header: "Sentence ID",
                                id: "textid",
                                accessor: d => d.textid,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['textid'] ];
                                    
                                    var retval = matchSorter(elems, filter.value);
        
                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    var sects = row.value.split(".").slice(-2)

                                    return <span>{sects.join(".")}</span>
                                }
                            },
                            {
                                Header: "Document Sections",
                                id: "sections",
                                accessor: d => d.sections,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = rows['sections'];
                                    
                                    var retval = matchSorter(elems, filter.value);
        
                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    return <span>{row.value.join(", ")}</span>
                                }
                            },
                            {
                                Header: "Ontology ID",
                                id: "obo_id",
                                accessor: d => d.obo_id,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['obo_id'] ];
                                    
                                    var retval = matchSorter(elems, filter.value);

                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    var allDisInfo = [];
                                
                                    var relInfo = [row.value];
    
                                    //console.log(relInfo);    
                                    for (var i = 0; i < relInfo.length; ++i)
                                    {   

                                        var linkID = relInfo[i].replace(":", "_");
                                        var linkID = linkID.replace('owl', '');

                                        if (linkID.indexOf('GO') !== -1)
                                        {
                                            allDisInfo.push(
                                                <span key={i} style={{display: "block"}}>                                        
                                                    <a href={"http://purl.obolibrary.org/obo/"+linkID}>{relInfo[i]}</a>
                                                </span>
                                                );
                                        } else if (linkID.indexOf('ATOL') !== -1)
                                        {
                                            allDisInfo.push(
                                                <span key={i} style={{display: "block"}}>                                        
                                                    <a href={"https://opendata.inra.fr/ATOL/page/"+linkID}>{relInfo[i]}</a>
                                                </span>
                                                );
                                        } else {
                                            allDisInfo.push(
                                                <span key={i} style={{display: "block"}}>                                        
                                                    {linkID}
                                                </span>
                                                );
                                        }
                                    }
    
                                    return <div>{allDisInfo}</div>;
                                }
                            },
                            {
                                Header: "Synonym",
                                id: "syn",
                                accessor: d => d.syn,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['syn'] ];
                                    
                                    var retval = matchSorter(elems, filter.value);

                                    return retval.length > 0;
                                },
                                Cell: (row) => {
                                    return <span>{row.value}</span>
                                }
                            },
                            {
                                Header: "Text",
                                id: "text",
                                accessor: d => d.text,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['text'] ];
                                    
                                    var retval = matchSorter(elems, filter.value);

                                    return retval.length > 0;
                                },
                                Cell: (row) => {
                                    return <span>{row.value}</span>
                                }
                            },
                            {
                                Header: "Context",
                                id: "context",
                                accessor: d => d.context,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['context'] ];
                                    
                                    var retval = matchSorter(elems, filter.value);

                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    var suffix = "...";
                                    if (row.value[row.value.length-1] == ".")
                                    {
                                        suffix = "";
                                    }

                                    return <span key={i} style={{display: "block"}}>...{row.value}{suffix}</span>
                                }
                            }
                        ]
                      }
                    ]}
                    defaultPageSize={10}
                    className="-striped -highlight"
                    
                  />);

                


            }


        }

        // PREPARE GRAPH HERE

        var graph= null;
        var interactionData = [];

        if (true)
        {
            var graphData = {
                nodes: [],
                links: []
            };

            var allAddedEdges = [];

            var allResPDFs = Object.keys(this.state.results);

            for (var pdfi = 0; pdfi < allResPDFs.length; ++pdfi)
            {
    
                var pdfData = this.state.results[allResPDFs[pdfi]];

                var sent2results = {}
    
                var pdfConcepts = Object.keys(pdfData);
    
                for (var j = 0; j < pdfConcepts.length; ++j)
                {
                    var contextData = pdfData[pdfConcepts[j]];

                    contextData.results.forEach(element => {

                        element["context_id"] = pdfConcepts[j];
                        
                        if (element.textid in sent2results)
                        {
                            sent2results[element.textid].push(element)
                        } else {
                            sent2results[element.textid] = [element]
                        }

                    });
                }

                var allSentIDs = Object.keys(sent2results);

                for (var i = 0; i < allSentIDs.length; ++i)
                {
                    
                    var context2res = {}
                    var sentid = allSentIDs[i];

                    sent2results[sentid].forEach(element => {

                        if (element["context_id"] in context2res)
                        {
                            context2res[element["context_id"]].push(element);
                        } else {
                            context2res[element["context_id"]] = [element];
                        }

                    });

                    var allFoundContexts = Object.keys(context2res);

                    if (allFoundContexts.length <= 1)
                    {
                        continue;
                    }

                    for (var ii=0; ii < allFoundContexts.length; ++ii)
                    {
                        for (var jj=ii+1; jj < allFoundContexts.length; ++jj)
                        {

                            var iContextRes = context2res[allFoundContexts[ii]];
                            var jContextRes = context2res[allFoundContexts[jj]];

                            
                            for (var iii=0; iii < iContextRes.length; ++iii)
                            {
                                for (var jjj=0; jjj < jContextRes.length; ++jjj)
                                {

                                    var resI = iContextRes[iii];
                                    var resJ = jContextRes[jjj];

                                    var foundSrcElems = graphData.nodes.filter((elem) => elem.id == resI.obo_id)
                                    if (foundSrcElems.length == 0)
                                    {

                                        var elemName = resI.syn + " (" + resI.obo_id + ")";

                                        if (resI.syn == resI.obo_id)
                                        {
                                            elemName = resI.syn;
                                        }
                                        // add node
                                        let nodeElem = {id: resI.obo_id, group: allFoundContexts[ii], name: elemName, idx: graphData.nodes.length, origname: resI.syn};
                                        graphData.nodes.push(nodeElem);
                                        foundSrcElems.push(nodeElem);
                                    }

                                    var foundTgtElems = graphData.nodes.filter((elem) => elem.id == resJ.obo_id)
                                    if (foundTgtElems.length == 0)
                                    {

                                        var elemName = resJ.syn + " (" + resJ.obo_id + ")";

                                        if (resJ.syn == resJ.obo_id)
                                        {
                                            elemName = resJ.syn;
                                        }

                                        // add node
                                        let nodeElem = {id: resJ.obo_id, group: allFoundContexts[jj], name: elemName, idx: graphData.nodes.length, origname: resJ.syn};
                                        graphData.nodes.push(nodeElem);
                                        foundTgtElems.push(nodeElem);
                                    }

                                    var sourceIdx = graphData.nodes.indexOf(foundSrcElems[0]);
                                    var targetIdx = graphData.nodes.indexOf(foundTgtElems[0]);
                                    
                                    var shortEdgeObj = {src: sourceIdx, tgt: targetIdx};
                                    if (allAddedEdges.indexOf(shortEdgeObj) >= 0)
                                    {
                                        continue;
                                    }


                                    if (((foundSrcElems[0].group == "Human Genes") && (foundTgtElems[0].group == "miRNAs")) || ((foundTgtElems[0].group == "Human Genes") && (foundSrcElems[0].group == "miRNAs")))
                                    {


                                        var newelem = {
                                            src: foundSrcElems[0],
                                            tgt: foundTgtElems[0],
                                            evidence_docs: [],
                                            sentence: [sentid]
                                        };

                                        if ((foundTgtElems[0].group == "Human Genes") && (foundSrcElems[0].group == "miRNAs"))
                                        {
                                            var tmp = newelem.src;
                                            newelem.src = newelem.tgt;
                                            newelem.tgt = tmp;
                                        }

                                        var foundInteractionData = interactionData.filter((elem) => elem.src.origname == newelem.src.origname && elem.tgt.origname == newelem.tgt.origname)

                                        if (foundInteractionData.length == 0)
                                        {
                                            interactionData.push(
                                                newelem
                                            )
                                        } else {
                                            foundInteractionData[0].sentence.push(sentid);
                                        }


                                    }

                                    allAddedEdges.push(shortEdgeObj);
                                    graphData.links.push({source: sourceIdx, target: targetIdx, group1: 1, group2: 1, group3: 1})

                                }
                            }



                        }
                    }


                }


                    allPDFTables.push(<h2 key={allPDFTables.length}>{this.resid2file[allResPDFs[pdfi]]} - Interaction Evidences</h2>)
                    allPDFTables.push(<ReactTable
                    key={allPDFTables.length}
                    data={this.state.interactionData}
                    filterable
                    defaultFilterMethod={(filter, row) =>
                        String(row[filter.id]) === filter.value}
                    columns={[
                        {
                        Header: "Found Concepts",
                        columns: [
                            {
                                Header: "Gene",
                                id: "srcid",
                                accessor: d => d.src.origname,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['srcid'] ];
                                    var retval = matchSorter(elems, filter.value);
        
                                    return retval.length > 0;
                                },
                                Cell: (row) => {   
                                    return <span>{row.value}</span>
                                }
                            },
                            {
                                Header: "miRNA",
                                id: "tgtid",
                                accessor: d => d.tgt.origname,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = [ rows['tgtid'] ];                                 
                                    var retval = matchSorter(elems, filter.value);
        
                                    return retval.length > 0;
                                },
                                Cell: (row) => {   
                                    return <span>{row.value}</span>
                                }
                            },
                            {
                                Header: "Sentences",
                                id: "sent_evidences",
                                accessor: d => d.sentence,
                                filterMethod: (filter, rows) => 
                                {                                   
                                    var allSects = [];

                                    rows['sentence'].value.forEach(element => {

                                        var sects = element.split(".").slice(-2)
                                        var sentSects = sects.join(".");
                                        allSects.push(sentSects);
                                    });
                                    
                                    var retval = matchSorter(allSects, filter.value);
        
                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    var allSects = [];

                                    row.value.sort().forEach(element => {

                                        var sects = element.split(".").slice(-2)
                                        var sentSects = sects.join(".");
                                        allSects.push(<span key={allSects.length}>{sentSects}</span>);
                                    });

                                    return <div style={{display: "grid"}}>{allSects}</div>
                                }
                            },
                            {
                                Header: "Evidences",
                                id: "doc_evidences",
                                accessor: d => d.evidence_docs,
                                filterMethod: (filter, rows) => 
                                {
                                    var elems = rows['evidence_docs'];
                                    
                                    var retval = matchSorter(elems, filter.value);
        
                                    return retval.length > 0;
                                },
                                Cell: (row) => {

                                    var allPmids = []

                                    row.value.sort().forEach(element => {
                                        allPmids.push(<a href={"https://www.ncbi.nlm.nih.gov/pubmed/"+element} target="_blank" key={allPmids.length}>{element}</a>)
                                    });
    
                                    return <div style={{display: "grid"}}>{allPmids}</div>
                                }
                            },
                        ]
                        }
                    ]}
                    defaultPageSize={10}
                    className="-striped -highlight"
                    
                    />);



            }

            console.log("Spongebob")
            console.log(graphData);
            console.log(this.state);

            if ((interactionData.length > 0) && (this.state.updateInteractions==true))
            {

                var self=this;
                var allpairs = [];

                interactionData.forEach(element => {
                    allpairs.push([element['src'].origname, element['tgt'].origname])
                });
                
                console.log("Attempting to fetch new details")
                console.log(allpairs)

                axios.post(config.getRestAddress() + "/details", {'detailsfor': allpairs}, config.axiosConfig)
                .then(function (response) {
        
                    console.log("Received TM Results")
                    console.log(response.data)

                    var newInteractData = [];
                    interactionData.forEach(element => {

                        var gene = element['src'].origname;
                        var mirna = element['tgt'].origname;

                        if ((gene in response.data) && (mirna in response.data[gene]))
                        {
                            element["evidence_docs"] = response.data[gene][mirna];
                        }

                        newInteractData.push(element);

                    })

                    console.log(newInteractData);
        
                    self.setState({interactionData: newInteractData, updateInteractions: false});
                })
                .catch(function (error) {
                  console.log(error);
                  self.setState({interactionData: [], updateInteractions: false});
                });

            }


            if (graphData.nodes.length > 0)
            {

                graph = <Card style={{marginBottom: "20px"}}>
                            <CardText>
                                <D3SVGParallelLinesGraph graph={graphData} graphInfo={{}}/>
                            </CardText>
                        </Card>;
            }

        }



        
        return (<div>
                    <Card style={{marginBottom: "20px"}}>
                        <CardText>
                            {allPDFTables}
                        </CardText>
                    </Card>
                    {graph}
                </div>
                
        );

    }

    /**
     * Render the component.
     */
    render() {

        var tmresult = this.makeTMResult();

        const dropzoneStyle = {
            width  : "100%",
            height : "100px",
            border : "1px dashed black"
        };

        return (

            <div>
            <Card style={{marginBottom: "20px"}}>
                <CardTitle
                    title="pdfMiRExplore - Welcome!"
                    subtitle="Check your literature for animal welfare concepts."
                />
                <CardText >
                    <p>Drop your PDFs and scan for miRNA-Target Interactions.</p>
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