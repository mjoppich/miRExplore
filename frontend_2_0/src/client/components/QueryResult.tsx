import axios from 'axios';
import {matchSorter} from 'match-sorter'
import * as React from 'react';

import ResultTable from "./ResultTable";

import config from '../config';
import EvidenceReportButton from './EvidenceReportButton';

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

 const orgTLC2Long = {
    'hsa': 'Homo sapiens',
    'mmu': 'Mus musculus'
    };

 export default class QueryResult extends React.Component<QueryResultProps, QueryResultState> {

    public static defaultProps: Partial<QueryResultProps> = {
        showFeatures: false,
        showInteractionGraph: false,
        showShortEvidenceTable: false,
        showEvidenceTable: true,
        searchQuery: null
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

                        var longOrg = orgTLC2Long[element] || element;
                        
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


            if ('modelanats' in this.props.docInfos)
            {
                rowData['modelanats_info'] = this.makeRelInfos(docIDs, this.props.docInfos['modelanats'])
            }

            data.push(rowData);   
        }

        var evidenceGraph = <div></div>;

        if (this.props.showInteractionGraph)
        {
            //evidenceGraph = <MEEvidenceGraph data={data}/>;
        }





        return (<div>
            <ResultTable data={data} />
        </div>)

    }

   

    /*

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

    */
}