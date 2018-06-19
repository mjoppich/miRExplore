import * as React from 'react'; 
import {Card, CardActions, CardHeader, CardText} from 'material-ui/Card';
import FlatButton from 'material-ui/FlatButton';
import CheckIcon from 'material-ui/svg-icons/action/done';
import DisagreeIcon from 'material-ui/svg-icons/action/thumb-down';


import Paper from 'material-ui/Paper';
import SelectedElements from '../components/SelectedElements';
import ACInput from '../components/AutoComplete';
import axios from 'axios';
import config from '../config';
import D3GraphViewer from '../components/D3GraphViewer';
import Toggle from 'material-ui/Toggle';
import OrganismChipAC from '../components/OrganismChipAC';
import EntityChipAC from '../components/EntityChipAC';
import LinearProgress from 'material-ui/LinearProgress';
import matchSorter from 'match-sorter';
import ReactTable from 'react-table';

import EvidenceReportButton from '../components/EvidenceReportButton';
import MEEvidenceGraph from '../components/MEEvidenceGraph';
import MEFeatureViewer from '../components/MEFeatureViewer';


import OboChipAC from '../components/OBOChipAC';

export interface QueryResultProps { 
    foundRelations: any,
    docInfos: any,
    searchWords: any,
    showGeneStructure: boolean,
    showInteractionGraph: boolean
};
export interface QueryResultState { 

 };
 class QueryResult extends React.Component<QueryResultProps, QueryResultState> {

    orgTLC2Long = {
        'hsa': 'Homo sapiens',
        'mmu': 'Mus musculus'
    };

    constructor(props) {
        super(props);

    }

    render ()
    {
        var self = this;
        console.log(this.props.foundRelations);
        console.log(this.props.docInfos);

        var data = [];
        
        var foundRels = this.props.foundRelations;

        for (var i = 0; i < foundRels.length; ++i)
        {
            var foundRel = foundRels[i];

            var rowData = {};
            rowData['rid'] = foundRel['rid'];
            rowData['lid'] = foundRel['lid'];

            rowData['ltype'] = foundRel['ltype'];
            rowData['rtype'] = foundRel['rtype'];

            rowData['evidences'] = foundRel['evidences'];

            var docIDs = [];
            var orgInfo = {};
            var allorgs = [];

            for (var e=0; e < foundRel['evidences'].length; ++e)
            {
                var ev = foundRel['evidences'][e];

                if ( docIDs.indexOf(ev['docid']) == -1 )
                {
                    docIDs.push(ev['docid']);
                }

                if ('orgs' in ev)
                {
                    ev['orgs'].forEach(element => {

                        console.log("org element");

                        var longOrg = this.orgTLC2Long[element] || element;
                        
                        if (allorgs.indexOf(longOrg) == -1)
                        {
                            allorgs.push(longOrg);
                        }

                        if (ev['docid'] in orgInfo)
                        {
                            if (orgInfo[ev['docid']].indexOf(element) == -1)
                            {
                                orgInfo[ev['docid']].push(element);
                            }
                        } else {
                            orgInfo[ev['docid']] = [element];
                        }
                    });
                }
            }

            rowData['docids'] = docIDs;
            rowData['orgs'] = orgInfo;
            rowData['allorgs'] = allorgs;


            if ('categories' in this.props.docInfos)
            {
                var catInfos = this.props.docInfos['categories'];
                var relInfo = {};

                for (var d=0; d < docIDs.length; ++d)
                {
                    var docid = docIDs[d];

                    var catInfos = this.props.docInfos['categories'];
                    if (!(docid in catInfos))
                    {
                        continue;
                    }

                    var docinfo = catInfos[docid];

                    for (var c=0; c < docinfo.length; ++c)
                    {
                        var di = docinfo[c];

                        if (di.termid in relInfo)
                        {
                            relInfo[di.termid].push( [docid, di] );
                        } else {
                            relInfo[di.termid] =  [[docid, di]];
                        }
                    }
                }

                rowData['categories_info'] = relInfo;
            }

            if ('messengers' in this.props.docInfos)
            {
                var catInfos = this.props.docInfos['messengers'];
                var relInfo = {};

                for (var d=0; d < docIDs.length; ++d)
                {
                    var docid = docIDs[d];

                    var catInfos = this.props.docInfos['messengers'];
                    if (!(docid in catInfos))
                    {
                        continue;
                    }

                    var docinfo = catInfos[docid];

                    for (var c=0; c < docinfo.length; ++c)
                    {
                        var di = docinfo[c];

                        if (di.termid in relInfo)
                        {
                            relInfo[di.termid].push( [docid, di] );
                        } else {
                            relInfo[di.termid] =  [[docid, di]];
                        }
                    }
                }

                rowData['messengers_info'] = relInfo;
            }

            
            data.push(rowData);   
        }

        var evidenceGraph = <div></div>;

        if (this.props.showInteractionGraph)
        {
            evidenceGraph = <MEEvidenceGraph data={data}/>;
        }

        
        return (
            <div>
              <ReactTable
                data={data}
                filterable
                defaultFilterMethod={(filter, row) =>
                  String(row[filter.id]) === filter.value}
                columns={[
                  {
                    Header: "Found Interaction",
                    columns: [
                      {
                        Header: "Cell",
                        id: "lid",
                        accessor: d => d.lid,
                        filterMethod: (filter, rows) => matchSorter(rows, filter.value, { keys: ["lid"] }),
                        Cell: (row) => {
                            return <span>{row.value}</span>
                        }
                      },
                      {
                        Header: "Cell",
                        id: "rid",
                        accessor: d => d.rid,
                        filterMethod: (filter, rows) => matchSorter(rows, filter.value, { keys: ["rid"] }),
                        filterAll: true
                      }
                    ]
                  },
                  {
                    Header: "Info",
                    columns: [
                      {
                        Header: "Document",
                        accessor: "docids",
                        Cell: (row) => {

                            var rows = [];
                            for (var i = 0; i < row.value.length; i++) {

                                if (i > 10)
                                {
                                    continue;
                                }
                                // note: we add a key prop here to allow react to uniquely identify each
                                // element in this array. see: https://reactjs.org/docs/lists-and-keys.html
                                rows.push(
                                    <span key={rows.length} style={{display: "block"}}>
                                         <a href={"https://www.ncbi.nlm.nih.gov/pubmed/"+row.value[i]}>{row.value[i]}</a>
                                    </span>
                                );
                            }

                            rows.push(
                                <span key={rows.length} style={{display: "block"}}>
                                    ... {row.value.length} Documents
                                </span>);

                            return <div>{rows}</div>;
                        }
                      },
                      {
                        Header: "Organisms",
                        accessor: "allorgs",
                        Cell: (row) => {

                            var rows = [];
                            for (var i = 0; i < row.value.length; i++) {
                                // note: we add a key prop here to allow react to uniquely identify each
                                // element in this array. see: https://reactjs.org/docs/lists-and-keys.html
                                rows.push(
                                    <span key={i} style={{display: "block"}}>
                                         {row.value[i]}
                                    </span>
                                );
                            }
                            return <div>{rows}</div>;
                        },
                        filterMethod: (filter, row) => {
                            var filterID = filter.id;
                            var rowData = row[filterID];

                            var retval = matchSorter(rowData, filter.value);
                            console.log(retval);

                            return retval.length > 0;
                        }
                      }
                    ]
                  },
                  {
                    Header: "Context",
                    columns: [
                      {
                        Header: "Categories",
                        id: "categories_info",
                        accessor: (d) => {
                            if ('categories_info' in d)
                            {
                                return d['categories_info'];
                            } else {
                                return null;
                            }
                        },
                        Cell: (row) => {

                            if (row.value == null)
                            {
                                return <div>N/A</div>
                            } else {

                                var allDisInfo = [];
                                
                                var relInfo = row.value;

                                console.log(relInfo);
                                var relKeys = Object.keys(relInfo);

                                for (var i = 0; i < relKeys.length; ++i)
                                {
                                    var termID = relKeys[i];
                                    var tinfos = relInfo[termID];

                                    var docIDs = [];
                                    
                                    for (var j=0; j < tinfos.length; ++j)
                                    {
                                        docIDs.push(tinfos[j][0]);
                                    }
                                    var tentry = tinfos[0][1];

                                    var linkID = tentry.termid.replace(":", "_");

                                    allDisInfo.push(
                                    <span key={i} style={{display: "block"}}>  
                                        {tentry.termname} ({docIDs.join(", ")})                                      
                                    </span>
                                    );
                                }

                                return <div>{allDisInfo}</div>;

                            }
                        },
                        filterMethod: (filter, row) => {
                            var filterID = filter.id;
                            var rowData = row[filterID];

                            var doids = Object.keys(rowData);

                            var allTerms = [];

                            for (var i=0; i < doids.length; ++i)
                            {
                                var doidentries = rowData[doids[i]];

                                for (var j=0; j < doidentries.length; ++j)
                                {
                                    var termname = doidentries[j][1]['termname'];

                                    if (allTerms.indexOf(termname) < 0)
                                    {
                                        allTerms.push(termname);
                                    }
                                }
                            }
                            
                            console.log("categories filter");
                            console.log(allTerms);
                            console.log(filter);

                            var retval = matchSorter(allTerms, filter.value);
                            console.log(retval);

                            return retval.length > 0;
                        }
                      },
                      {
                        Header: "Messenger",
                        id: "messengers_info",
                        accessor: (d) => {
                            if ('messengers_info' in d)
                            {
                                return d['messengers_info'];
                            } else {
                                return null;
                            }
                        },
                        Cell: (row) => {

                            if (row.value == null)
                            {
                                return <div>N/A</div>
                            } else {

                                var allDisInfo = [];
                                
                                var relInfo = row.value;

                                console.log(relInfo);
                                var relKeys = Object.keys(relInfo);

                                for (var i = 0; i < relKeys.length; ++i)
                                {
                                    var termID = relKeys[i];
                                    var tinfos = relInfo[termID];

                                    var docIDs = [];
                                    
                                    for (var j=0; j < tinfos.length; ++j)
                                    {
                                        docIDs.push(tinfos[j][0]);
                                    }
                                    var tentry = tinfos[0][1];

                                    var linkID = tentry.termid.replace(":", "_");

                                    allDisInfo.push(
                                    <span key={i} style={{display: "block"}}>                                        
                                        {tentry.termname} ({docIDs.join(", ")})
                                    </span>
                                    );
                                }

                                return <div>{allDisInfo}</div>;

                            }
                        },
                        filterMethod: (filter, row) => {
                            var filterID = filter.id;
                            var rowData = row[filterID];

                            var doids = Object.keys(rowData);

                            var allTerms = [];

                            for (var i=0; i < doids.length; ++i)
                            {
                                var doidentries = rowData[doids[i]];

                                for (var j=0; j < doidentries.length; ++j)
                                {
                                    var termname = doidentries[j][1]['termname'];

                                    if (allTerms.indexOf(termname) < 0)
                                    {
                                        allTerms.push(termname);
                                    }
                                }
                            }
                            
                            var retval = matchSorter(allTerms, filter.value);
                            return retval.length > 0;
                        }
                      }
                    ]
                  }
                ]}
                defaultPageSize={10}
                showPaginationTop
                showPaginationBottom
                className="-striped -highlight"
                SubComponent={row => {
                    console.log(row);

                    return (
                        <div>


                            
                            <div style={{ padding: "20px" }}>
                            <ReactTable
                            filterable
                            defaultFilterMethod={(filter, row) => String(row[filter.id]) === filter.value}
                                    data={row['original']['evidences']}
                                    columns={
                                        [

                                            {
                                                Header: "Evidence",
                                                columns: [
                                                  {
                                                    Header: "Document",
                                                    id:"documentinfo",
                                                    accessor: d => [d.docid, d.rel_sentence],
                                                    Cell: (row) => {
                                                        var docID = row.value[0];
                                                        var docSent = row.value[1];

                                                        var allHTMLElems = [];
                                                        allHTMLElems.push(<a key={allHTMLElems.length} href={"https://www.ncbi.nlm.nih.gov/pubmed/"+docID} target="_blank">{docID}</a>);
                                                        allHTMLElems.push(<span key={allHTMLElems.length}>{docSent}</span>);
                                                        
                                                        if ('orgs' in row['original'])
                                                        {
                                                            for (var i = 0; i < row['original']['orgs'].length; ++i)
                                                            {
                                                                var org = row['original']['orgs'][i];
                                                                var longOrg = this.orgTLC2Long[org] || org;

                                                                allHTMLElems.push(<span key={allHTMLElems.length}>{longOrg}</span>);

                                                            }
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

                                                        if ('orgs' in row)
                                                        {
                                                            for (var i = 0; i < row['orgs'].length; ++i)
                                                            {
                                                                var org = row['orgs'][i];
                                                                var longOrg = this.orgTLC2Long[org] || org;

                                                                allTerms.push(longOrg);
                                                            }
                                                        }
                                                        
                                                        var retval = matchSorter(allTerms, filter.value);
                                                        return retval.length > 0;
                                                    }
                                                  },
                                                  {
                                                    Header: "Relation",
                                                    id: "relation",
                                                    accessor: d => [d.lid, d.rid],
                                                    Cell: (row) => {
                                                        return <div style={{display: "grid"}}>
                                                            <span style={{whiteSpace: 'normal'}}>{row.value[0]}</span>
                                                            <span style={{whiteSpace: 'normal'}}>{row.value[1]}</span>
                                                        </div>;
                                                    }
                                                  },
                                                  {
                                                    Header: "Sentence",
                                                    id: "fsent",
                                                    accessor: d => d.sentence,
                                                    Cell: (row) => {

                                                        var tev = row['original'];

                                                        console.log(tev);

                                                        var tevHighlights = [];

                                                        if (tev['rel_pos'] != null)
                                                        {
                                                            tevHighlights.push([
                                                                tev['rel_pos'][0],
                                                                tev['rel_pos'][1],
                                                                "green"
                                                            ])
                                                        }
                                                
                                                        if (tev['lpos'] != null)
                                                        {
                                                            tevHighlights.push([
                                                                tev['lpos'][0],
                                                                tev['lpos'][1],
                                                                "blue"
                                                            ])
                                                        }
                                                
                                                        if (tev['rpos'] != null)
                                                        {
                                                            tevHighlights.push([
                                                                tev['rpos'][0],
                                                                tev['rpos'][1],
                                                                "red"
                                                            ])
                                                        }
                                                
                                                        if (tev['sentence'] != null)
                                                        {
                                                            // ;
                                                            return <div style={{whiteSpace: 'normal'}}>{this.highlightedText(tev['sentence'], tevHighlights)}</div>;
                                                        } else {
                                                            return <span></span>
                                                        }
                                                    }
                                                  }
                                                ]
                                              },
       
                                              {
                                                Header: "Info",
                                                columns: [{
                                                    Header: "Categories",
                                                    id: "categories_det_col",
                                                    accessor: (d) => {
                                                        return [d.docid, d.rel_sentence]
                                                    },
                                                    Cell: (row) => {

                                                        var docID = row.value[0];
                                                        var docSent = row.value[0];

                                                        console.log("Messenger Cell");
                                                        console.log(docID);
                                                        console.log(docSent);
                                                        console.log(this.props.docInfos);
                                                        console.log(row['original'])

                                                        if (!(docID in this.props.docInfos['categories']))
                                                        {
                                                            return <div>N/A</div>
                                                        }

                                                        var docCategories = this.props.docInfos['categories'][docID];                                                        

                                                        var allTerms = this.makeEvidenceTerms(docCategories, docID, docSent);
                            
                                                       return <div style={{display: 'grid'}}>
                                                            {allTerms}
                                                       </div>;
                                                    },
                                                    filterMethod: (filter, row) => {
                                                        var filterID = filter.id;
                                                        var rowData = row[filterID];

                            
                                                        var docID = rowData[0];
                                                        var sentID = rowData[1];

                                                        var docMessengers = this.props.docInfos['categories'][docID];       
                                                        
                                                        if (!docMessengers)
                                                        {
                                                            return false;
                                                        }

                                                        var allTerms = this.getPropTermsForDocID(docMessengers);
                                                        
                                                        var retval = matchSorter(allTerms, filter.value);
                                                        return retval.length > 0;
                                                    }
                                                
                                                },
                                                {
                                                    Header: "Messengers",
                                                    id: "messenger_det_col",
                                                    accessor: (d) => {
                                                        return [d.docid, d.rel_sentence]
                                                    },
                                                    Cell: (row) => {

                                                        var docID = row.value[0];
                                                        var docSent = row.value[0];

                                                        console.log("Messenger Cell");
                                                        console.log(docID);
                                                        console.log(docSent);
                                                        console.log(this.props.docInfos);
                                                        console.log(row['original'])

                                                        if (!(docID in this.props.docInfos['messengers']))
                                                        {
                                                            return <div>N/A</div>
                                                        }

                                                        var docMessengers = this.props.docInfos['messengers'][docID];                                                        

                                                        var allTerms = this.makeEvidenceTerms(docMessengers, docID, docSent);
                            
                                                       return <div style={{display: 'grid'}}>
                                                            {allTerms}
                                                       </div>;

                                                    },
                                                    filterMethod: (filter, row) => {
                                                        var filterID = filter.id;
                                                        var rowData = row[filterID];

                                                        console.log("Filter");
                                                        console.log(row);
                                                        console.log(rowData);
                            
                                                        var docID = rowData[0];
                                                        var sentID = rowData[1];

                                                        var docMessengers = this.props.docInfos['messengers'][docID];

                                                        if (!docMessengers)
                                                        {
                                                            return false;
                                                        }
                                                        
                                                        console.log(docMessengers);
                                                        var allTerms = this.getPropTermsForDocID(docMessengers);
                                                        
                                                        var retval = matchSorter(allTerms, filter.value);
                                                        return retval.length > 0;
                                                    }
                                                
                                                }]
                                              },
                                              {
                                                Header: "Feedback",
                                                columns: [{
                                                    Header: "Feedback",
                                                    id: "feedbackcol",
                                                    accessor: (d) => {
                                                        return d.docid
                                                    },
                                                    Cell: (row) => {
                            
                                                        if (row.value != null)
                                                        {
                                                            var tev = row['original'];
                                                            return <div style={{display: "grid"}}>
                                                                <EvidenceReportButton buttonDivStyle={{display: "grid"}} dataID={tev['data_id']} onAccept={() => self.reportEvidence(tev, true)} onDisagree={() => self.reportEvidence(tev, false)} />
                                                            </div>
                                                        } else {
                                                            return <div></div>;
                                                        }
                                                    }
                                                
                                                }]
                                              }
                                        ]
                                    }
                                    defaultPageSize={10}
                                    showPaginationTop
                                    showPaginationBottom
                                    className="-striped -highlight"

                                    />


                            </div>
                        </div>
                        )
                        }
                        }
              />
                {
                    evidenceGraph
                }
            </div>
          );

          /*
                                {self.prepareEvidences(row['original'])}



                                  <pre>{JSON.stringify(row, null, 2)}</pre>
                        <br />
                        <br />*/
    }

    getPropTermsForDocID(docMessengers)
    {
        var terms = [];

        for (var d=0; d < docMessengers.length; ++d)
        {
            var elem = docMessengers[d];

            var termName = elem['termname'];
            var termID = elem['termID'];

            if (terms.indexOf(termName) == -1)
            {
                terms.push(termName)
            }

            for (var i = 0; i < elem['evidences'].length; ++i)
            {
                var ev = elem['evidences'][i];
                var evSent = ev[0];

                var word = ev[3];
                var newelem = [termName, termID];
            }  
        }

        return terms;
    }

    makeEvidenceTerms(docMessengers, docID, docSent)
    {
        var occs = [];
                                                        var occsBySent = {};

                                                        for (var d=0; d < docMessengers.length; ++d)
                                                        {
                                                            var elem = docMessengers[d];
                                                            console.log(elem);

                                                            var termName = elem['termname'];
                                                            var termID = elem['termID'];

                                                            for (var i = 0; i < elem['evidences'].length; ++i)
                                                            {
                                                                var ev = elem['evidences'][i];
                                                                
                                                                console.log(ev[0]);
                                                                console.log(docSent);
                                                                
                                                                var evSent = ev[0];

                                                                var word = ev[3];
                                                                var newelem = [termName, termID];
    
                                                                var occsIndex = occs.indexOf(newelem);
                                                                if (occsIndex == -1)
                                                                {
                                                                    occsBySent[occs.length] = [[evSent, word]];
                                                                    occs.push( newelem );
                                                                } else {
                                                                    occsBySent[occsIndex].push([evSent, word]);
                                                                }


                                                            }
                                                        }
                                                        
                                                        var allTerms = []

                                                        for (var i = 0; i < occs.length; ++i)
                                                        {
                                                            var e = occs[i];
                                                            allTerms.push(<span key={2*i}>{e[0]} / {e[1]}</span>);
                                                            allTerms.push(<span key={2*i+1}>{occsBySent[i].join(", ")}</span>);
                                                        }

                                    return allTerms;
    }

    highlightedText( sentence, highlightAt)
    {

        var sortedHighlight = highlightAt.sort((pos1, pos2) => {
            if (pos1[0] < pos2[0])
            {
                return -1;
            }

            if (pos1[0] == pos2[0])
            {
                if (pos1[1] < pos2[1])
                {
                    return -1;
                } else {
                    return 1;
                }
            }

            return 1
        }).reverse();

        var allParts = [];
        var allText = "";

        var lastRight = sentence.length;
        var lastLeft = sentence.length;

        console.log(sortedHighlight);

        for (var i = 0; i < sortedHighlight.length; ++i)
        {
            var intStart = sortedHighlight[i][0];
            var intStop = sortedHighlight[i][1];
            var highColor = sortedHighlight[i][2];

            var rightPart = sentence.substr(intStop, lastLeft-intStop)
            var highlightPart = sentence.substr(intStart, intStop-intStart)

            allText = rightPart + allText

            allParts.push(<div key={allParts.length}>{rightPart}</div>);
            
            allParts.push(<span style={{color:highColor, fontWeight:"bold"}} key={allParts.length}>{highlightPart}</span>);

            lastRight = intStop;
            lastLeft = intStart;
        }

        allParts.push(<div key={allParts.length}>{sentence.substr(0, sortedHighlight[sortedHighlight.length-1][0])}</div>);

        allParts.reverse();

        return <span>{allParts}</span>
    }

    prepareEvidences( allInfo )
    {


        console.log(allInfo);
        var self = this;

        /*

        "docid": "25666935",
                    "gene_pos": [
                        64,
                        69
                    ],
                    "mirna_pos": [
                        19,
                        26
                    ],
                    "rel_category": "NEU",
                    "rel_direction": "MG",
                    "rel_direction_verb": "VMG",
                    "rel_negated": false,
                    "rel_pos": [
                        4,
                        15
                    ],
                    "rel_sentence": "25666935.2.6",
                    "rel_verb": "correlat",
                    "same_paragraph": true,
                    "same_sentence": true,
                    "sentence": "The correlation of miR-346 levels with the percentages of CD4(+)CXCR5(+)T cells and autoantibody levels were also analyzed."

            */

        var allEvs = allInfo['evidences'];
        var evStuff = [];

        var sortedEvs = allEvs.sort((lev, rev) => {
            if ((lev['data_source'] == 'pmid') && (rev['data_source'] != 'pmid'))
            {
                return -1;
            } else if ((lev['data_source'] == 'pmid') && (rev['data_source'] == 'pmid'))
            {
                return 0;
            } else  if ((lev['data_source'] != 'pmid') && (rev['data_source'] == 'pmid'))
            {
                return 1;
            }

            return 1;
        }).reverse();


        var beforeEvType = '';
        for (var i = 0; i < allEvs.length; ++i)
        {

            let tev = allEvs[i];

            //tev['lid'] = allInfo['lid']
            //tev['rid'] = allInfo['rid']

            if (tev['data_source'] == 'pmid')
            {

                if (beforeEvType != 'pmid')
                {
                    // add header
                    var headRow = <tr key={evStuff.length} style={{textAlign: 'left'}}>
                        <th>Found Relation</th>
                        <th>Verb-Model</th>
                        <th>Evidence Location</th>
                        <th>Organisms</th>
                        <th>DB Interaction</th>
                    </tr>;
            
                    evStuff.push(headRow);

                    beforeEvType = 'pmid';
                }

                var infoRows =this.prepareInfoRowPubmed(tev, evStuff.length, "");

                if (infoRows != null)
                {
                    evStuff = evStuff.concat(infoRows);
                }

            } else if (tev['data_source'] == 'mirecords')
            {
                if (beforeEvType != 'mirecords')
                {
                    // add header
                    var headRow = <tr key={evStuff.length}>
                    <th>Gene</th>
                    <th>miRNA</th>
                    <th>Data Source</th>
                    <th>Data Evidence</th>
                    <th>Report</th>
                    </tr>;
            
                    evStuff.push(headRow);

                    beforeEvType = 'mirecords';
                }


                var infoRows =this.prepareInfoRowMirecords(tev, evStuff.length);
                evStuff = evStuff.concat(infoRows);
            } else if(tev['data_source'] == 'miRTarBase')
            {

                if (beforeEvType != 'miRTarBase')
                {
                    // add header
                    var headRow = <tr key={evStuff.length}>
                    <th>Gene</th>
                    <th>miRNA</th>
                    <th>Data Source</th>
                    <th>Data Evidence</th>
                    <th>Report</th>
                    </tr>;
            
                    evStuff.push(headRow);

                    beforeEvType = 'miRTarBase';
                }

                var infoRows =this.prepareInfoRowMirTarBase(tev, evStuff.length);
                evStuff = evStuff.concat(infoRows);
            }



        }

        return <table style={{ tableLayout: "auto", width: "100%"}}>
                <tbody>
                    {evStuff}
                </tbody>
            </table>;

    }

    prepareInfoRowMirTarBase(tev, idx)
    {

        /*

                {
                    "data_id": "MIRT054715",
                    "data_source": "miRTarBase",
                    "docid": "24141785",
                    "exp_support": [
                        "Luciferase reporter assay",
                        "Western blot"
                    ],
                    "functional_type": "Functional MTI",
                    "lid": "CXCR4",
                    "ltype": "gene",
                    "organism": "Homo sapiens",
                    "rid": "miR-miR-9-5p",
                    "rtype": "mirna"
                },

                */

        var outRows = [];

        var infoRow = <tr key={idx+1}>
                        <td>{tev['lid']} ({tev['ltype']})</td>
                        <td>{tev['rid']} ({tev['rtype']})</td>
                        <td><a href={"http://mirtarbase.mbc.nctu.edu.tw/php/detail.php?mirtid="+tev['data_id']}>{tev['data_id']} ({tev['data_source']})</a></td>
                        <td>
                            <span style={{display:"block"}}><a href={"https://www.ncbi.nlm.nih.gov/pubmed/"+tev['docid']}>{tev['docid']}</a></span>
                            <span style={{display:"block"}}>{tev['functional_type']}</span>
                        </td>
                        <td>
                            {
                                tev['exp_support'].map((exptype, nsi) => <span key={nsi}  style={{display: "block"}}>{exptype}</span>)
                            }
                        </td>
                      </tr>;

        outRows.push(infoRow);

        return outRows;
    }

    prepareInfoRowMirecords(tev, idx)
    {

        /*

                        {
                    "data_source": "mirecords",
                    "docid": "18568019",
                    "lid": "CXCR4",
                    "ltype": "gene",
                    "rid": "hsa-miR-146a",
                    "rtype": "mirna"
                }

                */

        var outRows = [];

        var self=this;



        var infoRow = <tr key={idx}>
                        <td>{tev['lid']}, {tev['ltype']}</td>
                        <td>{tev['rid']}, {tev['rtype']}</td>
                        <td>{tev['data_id']} ({tev['data_source']})</td>
                        <td>{tev['docid']}</td>
                        <td>
                        <FlatButton label="Accept Entry" onClick={() => self.reportEvidence(tev, true)} icon={<CheckIcon/>}/>
                        <FlatButton label="Disagree Entry" backgroundColor="secondary" onClick={() => self.reportEvidence(tev, false)} icon={<DisagreeIcon/>}/>
                        </td>
                      </tr>;

        outRows.push(infoRow);

        return outRows;
    }

    prepareInfoRowPubmed(tev, idx, outlinkBase)
    {

        var outRows = [];

        var self=this;

        var relDirection;
        if (tev['rel_direction_verb'] != null)
        {
            relDirection = tev['rel_direction_verb'];
        } else {
            relDirection = tev['rel_direction'];
        }

        var relNegated = "";
        if (tev['rel_negated'] == true)
        {
            relNegated = ", negated";
        }

        var relLocation = ""
        
        if (tev['rel_sentence'] != null)
        {
            relLocation = tev['rel_sentence'] + " (same sentence)";
        } else {
            if (tev['same_paragraph'])
            {
                relLocation = "same paragraph";
            }
        }

        let acceptColor = ""; // -> #D0E2BF
        let disagreeColor = ""; // -> #e2bfbf

        var orgInfo = "";
        
        if ('orgs' in tev)
        {
            orgInfo = tev['orgs'].map((d) => self.orgTLC2Long[d] || d).join(", ");
        }

        var infoRow = <tr key={idx}>
                        <td>{tev['rel_verb']}, {tev['rel_category']}</td>
                        <td>{relDirection + relNegated}</td>
                        <td>{relLocation}</td>
                        <td>{orgInfo}</td>
                        <td rowSpan={2}>

                        <EvidenceReportButton buttonDivStyle={{}} dataID={tev['data_id']} onAccept={() => self.reportEvidence(tev, true)} onDisagree={() => self.reportEvidence(tev, false)} />

                        </td>
                        </tr>;

        outRows.push(infoRow);

        var tevHighlights = [];

        if (tev['rel_pos'] != null)
        {
            tevHighlights.push([
                tev['rel_pos'][0],
                tev['rel_pos'][1],
                "green"
            ])
        }

        if (tev['lpos'] != null)
        {
            tevHighlights.push([
                tev['lpos'][0],
                tev['lpos'][1],
                "blue"
            ])
        }

        if (tev['rpos'] != null)
        {
            tevHighlights.push([
                tev['rpos'][0],
                tev['rpos'][1],
                "red"
            ])
        }


        var sentRow = null;
        if (tev['sentence'] != null)
        {
            sentRow = <tr key={idx+1}><td colSpan={4}>{this.highlightedText(tev['sentence'], tevHighlights)}</td></tr>
        } else {
            return null;
            //sentRow = <tr key={idx+1}><td colSpan={4}>{tev['sentence']}</td></tr>
        }

        outRows.push(sentRow);

        return outRows;

}

    makeSearchTerms(searchWords)
    {
        var searchData = {};
        for (var i = 0; i <searchWords.length; ++i)
        {
            let elem = searchWords[i];

            let elemGroup = elem['group'];
            let elemName = elem['name'];
            
            if (elemGroup in searchData)
            {
                searchData[elemGroup].push(elemName);
            } else {
                searchData[elemGroup] = [elemName];
            }

        }

        return searchData;
            
    }

    reportEvidence(ev, accept)
    {
        console.log("Evidence Reported");
        console.log(accept);
        console.log(ev);

        var sendData = ev;
        sendData['approve'] = accept;
        sendData['search_terms'] = this.makeSearchTerms(this.props.searchWords)

        if ('lontid' in sendData)
        {
            var oldontid = sendData['lontid'];
            sendData['londit'] = sendData['lid']
            sendData['lid'] = oldontid;
        }

        if ('rontid' in sendData)
        {
            var oldontid = sendData['rontid'];
            sendData['rondit'] = sendData['rid']
            sendData['rid'] = oldontid;
        }

        axios.post(config.getRestAddress() + "/relation_feedback",sendData, config.axiosConfig)
        .then(function (response) {
          console.log(response.data)

        })
        .catch(function (error) {
          console.log(error);
        });
    }
 }
  
export interface QueryComponentProps { key: number};
export interface QueryComponentState { 
    selectedElements: Array<any>,
    selectedOrganisms: Array<any>,
    selectedCategories: Array<any>,
    selectedMessengers: Array<any>,
    interactions: any,
    showGeneStructure: boolean,
    showInteractionGraph:boolean,
    predictiveInteractions: boolean,
    queryStarted: boolean
 };

class QueryComponent extends React.Component<QueryComponentProps, QueryComponentState> {
    constructor(props) {
        super(props);

    }

    componentWillMount()
    {
        this.setState({selectedElements: [], predictiveInteractions:true,showInteractionGraph: false, showGeneStructure: false, selectedOrganisms: [], selectedMessengers: [], selectedCategories: []});
    }

    newElementSelected( newElement )
    {
        this.state.selectedElements.push(newElement);
        this.setState({selectedElements: this.state.selectedElements});
    }

    elementClicked( elementText )
    {
        console.log("Clicked on: " + elementText)
    }

    deleteElement( elementText, i )
    {
        var idx = this.state.selectedElements.indexOf(elementText);

        if (idx >= 0)
        {
            this.state.selectedElements.splice(idx, 1);
        }

        this.setState({selectedElements: this.state.selectedElements})
    }

    prepareResults()
    {
        var self = this;

        let sendData = {};

        console.log("Selected Elements in Explore")
        console.log(this.state.selectedElements);
        console.log(this.state);
        console.log(sendData);

        if (!(this.state.selectedElements))
        {
            return;
        }

        for (var i = 0; i < this.state.selectedElements.length; ++i)
        {
            sendData['elements'] = this.state.selectedElements;   
        }

        if ((this.state.selectedMessengers) && (this.state.selectedMessengers.length > 0))
        {
            sendData['messengers'] = this.state.selectedMessengers;
        }

        if ((this.state.selectedCategories) && (this.state.selectedCategories.length > 0))
        {
            sendData['categories'] = this.state.selectedCategories;
        }

        if ((this.state.selectedOrganisms)&&(this.state.selectedOrganisms.length > 0))
        {
            sendData['organisms'] = this.state.selectedOrganisms;
        }

        console.log(sendData);

        this.setState({queryStarted: true});

        axios.post(config.getRestAddress() + "/find_interactions",sendData, config.axiosConfig)
        .then(function (response) {
          console.log(response.data)

          self.setState({interactions: response.data, queryStarted: false})

        })
        .catch(function (error) {
          console.log(error);
          self.setState({interactions: {}, queryStarted: false});
        });
    }



    render()
    {
        
/*


 <Toggle
                label="Include Predictive Interactions"
                defaultToggled={true}
                toggled={this.state.predictiveInteractions}
                onToggle={(event, newValue) => this.setState({predictiveInteractions: !this.state.predictiveInteractions})}
                />
                <Toggle
                label="Show Gene Structure"
                defaultToggled={false}
                toggled={this.state.showGeneStructure}
                onToggle={(event, newValue) => this.setState({showGeneStructure: !this.state.showGeneStructure})}
                />

                */


        var alignResults = [];

        if (this.state.queryStarted)
        {
            alignResults.push(<LinearProgress key={alignResults.length}/>);
            alignResults.push(<p key={alignResults.length}>Retrieving Results</p>);

        }

        if ((this.state.interactions == null) || (this.state.interactions.length == 0))
        {

            if (!this.state.queryStarted)
            {
                alignResults.push(<p key={0}>No Result Available for your query.</p>)   ;
            }
            //alignResults.push(<pre key={1}>{JSON.stringify(this.state.selectedOrganisms, null, 2)}</pre>)   ;

        } else {

            alignResults.push(<QueryResult key={0} showInteractionGraph={this.state.showInteractionGraph} showGeneStructure={this.state.showGeneStructure} searchWords={this.state.selectedElements} foundRelations={this.state.interactions["rels"]} docInfos={this.state.interactions["pmidinfo"]}/>)

            //alignResults.push(<pre key={0}>{JSON.stringify(this.state.interactions, null, 2)}</pre>)   ;
            //alignResults.push(<pre key={1}>{JSON.stringify(this.state.selectedOrganisms, null, 2)}</pre>)   ;
            //var alignKeys = Object.keys(this.state.alignments);   
            
            //<EntityChipAC onValueChange={(newvalues) => {console.log("onVC called"); this.setState({selectedElements: newvalues})}} />
        }

        return (<Card style={{marginBottom: "20px"}}>
                <CardHeader
                title="Search Homology Entries"
                subtitle="Search by Gene/Protein ID"
                />
                <CardText>

                    <div>
                        

                        <OboChipAC
                            url="autocomplete"
                            floatText="Cell-Type"
                            hintText="Enter Cell-Type"
                            onValueChange={(newvalues) => this.setState({selectedElements: newvalues})
                        }/>

                        <OboChipAC
                            url="organisms"
                            floatText="Organism"
                            hintText="Enter organism name"
                            onValueChange={(newvalues) => this.setState({selectedOrganisms: newvalues})
                        }/>


                        <OboChipAC
                            url="category_ac"
                            floatText="Categories"
                            hintText="Enter category name"
                            onValueChange={(newvalues) => this.setState({selectedCategories: newvalues})
                        }/>

                        <OboChipAC
                            url="messengers_ac"
                            floatText="Messengers"
                            hintText="Enter Messenger-Term here"
                            onValueChange={(newvalues) => this.setState({selectedMessengers: newvalues})
                        }/>

               
                <Toggle
                label="Show Interaction Graph"
                defaultToggled={false}
                toggled={this.state.showInteractionGraph}
                onToggle={(event, newValue) => this.setState({showInteractionGraph: !this.state.showInteractionGraph})}
                />
                        <FlatButton label="Query specified Elements" onClick={() => this.prepareResults()}/>
                    </div>

                    <div>
                        {
                            alignResults
                        }
                        
                    </div>

                </CardText>
            </Card>);
    }
};

export interface ExplorePageProps { };
export interface ExplorePageState { queriesStored: number};

export class ExploreMainPage extends React.Component<ExplorePageProps, ExplorePageState> {

    allQueries = [];
    /**
     * Class constructor.
     */
    constructor(props) {
        super(props);
    }

    /**
     * This method will be executed after initial rendering.
     */
    componentWillMount() {

        this.setState({queriesStored: 0});

    }

    newQuery()
    {
        this.allQueries.push(<QueryComponent key={this.allQueries.length}/>);
        this.setState({queriesStored: this.allQueries.length});
    }

    clearQueries()
    {
        this.allQueries = [];
        this.setState({queriesStored: this.allQueries.length});
    }

    /**
     * Render the component.
     */
    render() {


        return (

            <div>

                <Card style={{marginBottom: "20px"}}>
                    <CardHeader
                    title="Create Query"
                    subtitle="Manage your queries"
                    actAsExpander={true}
                    showExpandableButton={true}
                    />
                    <CardActions>
                    <FlatButton label="New Query" onClick={this.newQuery.bind(this)}/>
                    <FlatButton label="Clear" onClick={this.clearQueries.bind(this)}/>
                    </CardActions>
                    <CardText expandable={true}>

                        <p>Some explanation on how to use this view!</p>

                    </CardText>
                </Card>
                
                <div>
                {this.allQueries}
                </div>

            </div>

        );
    }

}
