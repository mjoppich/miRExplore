import axios from 'axios';
import matchSorter from 'match-sorter';
import * as React from 'react';
import ReactTable from 'react-table';

import EvidenceReportButton from '../components/EvidenceReportButton';
import MEEvidenceGraph from '../components/MEEvidenceGraph';
import config from '../config';
import NISankeyChart from './NISankeyChart';
import NIInteractionNetwork from './NIInteractionNetwork';


export interface QueryResultProps { 
    foundRelations: any,
    docInfos: any,
    searchWords: any,
    showSankeyChart: boolean,
    showInteractionGraph: boolean,
    showEvidenceTable: boolean,
    showShortEvidenceTable: boolean,
    searchQuery: any,
    obolevels?: any
};
export interface QueryResultState { 

 };

 export default class QueryResult extends React.Component<QueryResultProps, QueryResultState> {

    public static defaultProps: Partial<QueryResultProps> = {
        showSankeyChart: false,
        showInteractionGraph: false,
        showShortEvidenceTable: false,
        showEvidenceTable: true,
        searchQuery: null,
        obolevels: {cells: 3, messengers: 3}
        };

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

    var evidenceGraph = [];

    //             evidenceGraph.push( <MEEvidenceGraph key={evidenceGraph.length} data={data}/>);


        if (this.props.showInteractionGraph)
        {
            evidenceGraph.push( <NIInteractionNetwork key={evidenceGraph.length} obolevels={this.props.obolevels} data={this.props.searchQuery} graphtitle={"Cell-Cell Interactions"}/>);

        }

        if (this.props.showSankeyChart)
        {
            evidenceGraph.push( <NISankeyChart key={evidenceGraph.length} obolevels={this.props.obolevels} data={this.props.searchQuery} graphtitle={"Cell-Cell Interactions to Category and Location"}/>);
        }

        
        return (
            <div>
                {evidenceGraph}
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
                            return this.makeSubcomponent(row);
                        }
                        }
              />

            </div>
          );

    }
    

    makeSubcomponent(row)
    {
        var self=this;
        return (<div>


                            
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
                                                        allTerms.push(rowData[0]);
                                                    
                                                        
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
                                                    },
                                                    filterMethod: (filter, row) => {
                                                        var filterID = filter.id;
                                                        var rowData = row[filterID];

                                                        var allTerms = [];

                                                        allTerms.push(rowData[0]);
                                                        allTerms.push(rowData[1]);

                                                        var retval = matchSorter(allTerms, filter.value);

                                                        console.log("rel Filter")
                                                        console.log(retval)
                                                        console.log(filter.value)
                                                        console.log(allTerms)
                                                        
                                                        return retval.length > 0;
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
                                                    },
                                                    filterMethod: (filter, row) => {
                                                        var filterID = filter.id;
                                                        var rowData = row[filterID];

                                                        var allTerms = rowData.split(" ");
                                                        
                                                        var retval = matchSorter(allTerms, filter.value);

                                                        console.log("Sentence Filter")
                                                        console.log(allTerms)
                                                        console.log(retval)
                                                        console.log(filter.value)

                                                        return retval.length > 0;
                                                    }
                                                  }
                                                ]
                                              },
                                              {
                                                Header: "Messengers",
                                                id: "messenger_det_col",
                                                accessor: (d:any) => {
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
                                                
                                                }
                                                ]
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
                        );
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