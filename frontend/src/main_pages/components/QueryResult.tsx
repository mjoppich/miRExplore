import axios from 'axios';
import matchSorter from 'match-sorter';
import * as React from 'react';
import ReactTable from 'react-table';

import EvidenceReportButton from '../components/EvidenceReportButton';
import MEEvidenceGraph from '../components/MEEvidenceGraph';
import MEFeatureViewer from '../components/MEFeatureViewer';

import config from '../config';



export interface QueryResultProps { 
    foundRelations: any,
    docInfos: any,
    searchWords: any,
    showFeatures: boolean,
    showInteractionGraph: boolean,
    showEvidenceTable: boolean,
    showShortEvidenceTable: boolean,
    searchQuery: any,
};
export interface QueryResultState { 

 };

 export default class QueryResult extends React.Component<QueryResultProps, QueryResultState> {

    public static defaultProps: Partial<QueryResultProps> = {
        showFeatures: false,
        showInteractionGraph: false,
        showShortEvidenceTable: false,
        showEvidenceTable: true,
        searchQuery: null
        };

    
    orgTLC2Long = {
        'hsa': 'Homo sapiens',
        'mmu': 'Mus musculus'
    };

    constructor(props) {
        super(props);

    }

    makeRelInfos(docIDs, docInfos)
    {
        var relInfo = {};

        for (var d=0; d < docIDs.length; ++d)
        {
            var docid = docIDs[d];

            if (!(docid in docInfos))
            {
                continue;
            }

            var docinfo = docInfos[docid];

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

        return relInfo;
    }


    render ()
    {
        var self = this;
        //console.log(this.props.foundRelations);
        //console.log(this.props.docInfos);

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

                        //console.log("org element");

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


            if ('disease' in this.props.docInfos)
            {
                var catInfos = this.props.docInfos['disease'];
                var relInfo = {};

                for (var d=0; d < docIDs.length; ++d)
                {
                    var docid = docIDs[d];

                    var catInfos = this.props.docInfos['disease'];
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

                rowData['disease_info'] = relInfo;
            }

            if ('go' in this.props.docInfos)
            {
                var catInfos = this.props.docInfos['go'];
                var relInfo = {};

                for (var d=0; d < docIDs.length; ++d)
                {
                    var docid = docIDs[d];

                    var catInfos = this.props.docInfos['go'];
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

                rowData['go_info'] = relInfo;
            }

            if ('cells' in this.props.docInfos)
            {
                rowData['cells_info'] = this.makeRelInfos(docIDs, this.props.docInfos['cells'])
            }

            if ('ncits' in this.props.docInfos)
            {
                rowData['ncits_info'] = this.makeRelInfos(docIDs, this.props.docInfos['ncits'])
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
                        Header: "Gene/lncRNA",
                        id: "lid",
                        accessor: d => d.lid,
                        filterMethod: (filter, rows) => matchSorter(rows, filter.value, { keys: ["lid"] }),
                        Cell: (row) => {

                            //https://www.genecards.org/cgi-bin/carddisp.pl?gene=RUNX3
                            if (row.original.ltype == "gene")
                            {
                               return (<span style={{display: "block"}}>
                               <a href={"https://www.genecards.org/cgi-bin/carddisp.pl?gene="+row.value}>{row.value}</a>
                            </span>); 
                            }

                            return <span>{row.value}</span>
                        }
                      },
                      {
                        Header: "miRNA/lncRNA",
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
                                // note: we add a key prop here to allow react to uniquely identify each
                                // element in this array. see: https://reactjs.org/docs/lists-and-keys.html
                                rows.push(
                                    <span key={i} style={{display: "block"}}>
                                         <a href={"https://www.ncbi.nlm.nih.gov/pubmed/"+row.value[i]}>{row.value[i]}</a>
                                    </span>
                                );
                            }
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
                        Header: "Disease",
                        id: "disease_info",
                        accessor: (d) => {
                            if ('disease_info' in d)
                            {
                                return d['disease_info'];
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

                                //console.log(relInfo);
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
                                        <a href={"http://purl.obolibrary.org/obo/"+linkID}>{tentry.termname} ({docIDs.join(", ")})</a>
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
                            
                            //console.log("disease filter");
                            //console.log(allTerms);
                            //console.log(filter);

                            var retval = matchSorter(allTerms, filter.value);
                            //console.log(retval);

                            return retval.length > 0;
                        }
                      },
                      {
                        Header: "Protein Class",
                        id: "ncit_info",
                        accessor: (d) => {
                            if ('ncits_info' in d)
                            {
                                return d['ncits_info'];
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

                                //console.log(relInfo);
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
                                        <a href={"http://purl.obolibrary.org/obo/"+linkID}>{tentry.termname} ({docIDs.join(", ")})</a>
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
                            
                            //console.log("disease filter");
                            //console.log(allTerms);
                            //console.log(filter);

                            var retval = matchSorter(allTerms, filter.value);
                            //console.log(retval);

                            return retval.length > 0;
                        }
                      },
                      {
                        Header: "Gene Ontology",
                        id: "go_info",
                        accessor: (d) => {
                            if ('go_info' in d)
                            {
                                return d['go_info'];
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

                                //console.log(relInfo);
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
                                        <a href={"http://purl.obolibrary.org/obo/"+linkID}>{tentry.termname} ({docIDs.join(", ")})</a>
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
                            
                            //console.log("go filter");
                            //console.log(allTerms);
                            //console.log(filter);

                            var retval = matchSorter(allTerms, filter.value);
                            //console.log(retval);

                            return retval.length > 0;
                        }
                      },
                      {
                        Header: "Cellline/FMA",
                        id: "cells_info",
                        accessor: (d) => {
                            if ('cells_info' in d)
                            {
                                return d['cells_info'];
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

                                //console.log(relInfo);
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
                                        <a href={"http://purl.obolibrary.org/obo/"+linkID}>{tentry.termname} ({docIDs.join(", ")})</a>
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
                            
                            //console.log("cell filter");
                            //console.log(allTerms);
                            //console.log(filter);

                            var retval = matchSorter(allTerms, filter.value);
                            //console.log(retval);

                            return retval.length > 0;
                        }
                      }
                    ]
                  }
                ]}
                defaultPageSize={10}
                className="-striped -highlight"
                SubComponent={row => {
                    //console.log(row);

                console.log("subcomponent");
                console.log(self.props);
                console.log(row);

                var geneStructure = [<span key={0}></span>];

                    if (self.props.showFeatures)
                    {
                        geneStructure = [];

                        var allGeneLncs = [];

                        if (row.original['ltype'] == 'gene')
                        {
                            allGeneLncs.push(row.original['lid'])
                        } else if (row.original['rtype'] == 'gene') {
                            allGeneLncs.push(row.original['rid'])
                        }

                        console.log("Creating Feature Viewer for")
                        console.log(allGeneLncs)

                        for (var i = 0; i < allGeneLncs.length; ++i)
                        {
                            let geneID = allGeneLncs[i];

                            geneStructure.push(<MEFeatureViewer key={geneStructure.length} featureID={geneID}/>)
                        }

                    }

                    return (
                        <div>
                            <div>
                                {geneStructure}
                            </div>
                            <div style={{ padding: "20px" }}>

                                {self.prepareEvidences(row['original'])}

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
                                  <pre>{JSON.stringify(row, null, 2)}</pre>
                        <br />
                        <br />*/
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

        var lastRight = sentence.length;
        var lastLeft = sentence.length;

        //console.log(sortedHighlight);

        for (var i = 0; i < sortedHighlight.length; ++i)
        {
            var intStart = sortedHighlight[i][0];
            var intStop = sortedHighlight[i][1];
            var highColor = sortedHighlight[i][2];

            var rightPart = sentence.substr(intStop, lastLeft-intStop)
            var highlightPart = sentence.substr(intStart, intStop-intStart)

            allParts.push(<span key={allParts.length}>{rightPart}</span>);
            allParts.push(<span style={{color:highColor, fontWeight:"bold"}} key={allParts.length}>{highlightPart}</span>);

            lastRight = intStop;
            lastLeft = intStart;
        }

        allParts.push(<span key={allParts.length}>{sentence.substr(0, sortedHighlight[sortedHighlight.length-1][0])}</span>);

        allParts.reverse();

        return <span>{allParts}</span>
    }

    prepareEvidences( allInfo )
    {


        //console.log(allInfo);
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
                evStuff = evStuff.concat(infoRows);
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

                        <EvidenceReportButton dataID={tev['data_id']} onAccept={() => self.reportEvidence(tev, true)} onDisagree={() => self.reportEvidence(tev, false)} />

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

                        <EvidenceReportButton dataID={tev['data_id']} onAccept={() => self.reportEvidence(tev, true)} onDisagree={() => self.reportEvidence(tev, false)} />

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

        var sentRow = <tr key={idx+1}><td colSpan={4}>{""}</td></tr>;
        
        if (tev.sentence)
        {
            sentRow = <tr key={idx+1}><td colSpan={4}>{this.highlightedText(tev['sentence'], tevHighlights)}</td></tr>;
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

        axios.post(config.getRestAddress() + "/relation_feedback",sendData, config.axiosConfig)
        .then(function (response) {
          console.log(response.data)

        })
        .catch(function (error) {
          console.log(error);
        });
    }
}