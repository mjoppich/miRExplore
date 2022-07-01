
import * as React from 'react';
import Box from '@mui/material/Box';
import Collapse from '@mui/material/Collapse';
import IconButton from '@mui/material/IconButton';
import Table from '@mui/material/Table';
import TableBody from '@mui/material/TableBody';
import TableCell from '@mui/material/TableCell';
import TableContainer from '@mui/material/TableContainer';
import TableHead from '@mui/material/TableHead';
import TableRow from '@mui/material/TableRow';
import Typography from '@mui/material/Typography';
import Paper from '@mui/material/Paper';
import KeyboardArrowDownIcon from '@mui/icons-material/KeyboardArrowDown';
import KeyboardArrowUpIcon from '@mui/icons-material/KeyboardArrowUp';

import EvidenceReportButton from './EvidenceReportButton';

const orgTLC2Long = {
    'hsa': 'Homo sapiens',
    'mmu': 'Mus musculus'
    };

function foldOntologyData(indata:any, datakey:any, linkPrefix:any)
{
    if (! (datakey in indata))
    {
        return ""
    }

    var data = indata[datakey]

    var allDisInfo = [];                    

    //console.log(relInfo);
    var relKeys = Object.keys(data);

    var maxCount = 20;
    var iterateCount = Math.min(relKeys.length, maxCount);


    for (var i = 0; i < iterateCount; ++i)
    {
        var termID = relKeys[i];
        var tinfos = data[termID];

        var docIDs = [];
        
        for (var j=0; j < tinfos.length; ++j)
        {
            docIDs.push(tinfos[j][0]);
        }
        var tentry = tinfos[0][1];

        var linkID = tentry.termid.replace(":", "_");

        allDisInfo.push(
        <span key={i} style={{display: "block"}}>                                        
            <a href={linkPrefix+linkID}>{tentry.termname} ({docIDs.join(", ")})</a>
        </span>
        );
    }

    if (iterateCount != relKeys.length)
    {
        allDisInfo.push(<span key={allDisInfo.length} style={{display: "block"}}>{relKeys.length-iterateCount} elements not shown.</span>)
    }

    return <div>{allDisInfo}</div>;
}

function foldListData(data:any, linkPrefix?:any)
{

    if (data == null)
    {
        return ""
    }

    var rows = [];
    for (var i = 0; i < data.length; i++) {
        // note: we add a key prop here to allow react to uniquely identify each
        // element in this array. see: https://reactjs.org/docs/lists-and-keys.html

        if (linkPrefix)
        {
            rows.push(
                <span key={i} style={{display: "block"}}>
                     <a href={"https://www.ncbi.nlm.nih.gov/pubmed/"+data[i]}>{data[i]}</a>
                </span>
            );
        } else {
            rows.push(
                <span key={i} style={{display: "block"}}>
                    {data[i]}
                </span>
            );
        }
        
    }
    return <div>{rows}</div>;

}


function Row(props: { row: any }) {
    const { row } = props;
    const [open, setOpen] = React.useState(false);
  
    return (
      <React.Fragment>
        <TableRow sx={{ '& > *': { borderBottom: 'unset' } }}>
          <TableCell style={{ verticalAlign: 'top' }}>
            <IconButton
              aria-label="expand row"
              size="small"
              onClick={() => setOpen(!open)}
            >
              {open ? <KeyboardArrowUpIcon /> : <KeyboardArrowDownIcon />}
            </IconButton>
          </TableCell>
          <TableCell component="th" scope="row" style={{ verticalAlign: 'top' }}>{row.lid}</TableCell>
          <TableCell align="right" style={{ verticalAlign: 'top' }}>{row.rid}</TableCell>
          <TableCell align="right" style={{ verticalAlign: 'top' }}>{foldListData(row.docids, "https://www.ncbi.nlm.nih.gov/pubmed/")}</TableCell>
          <TableCell align="right" style={{ verticalAlign: 'top' }}>{foldListData(row.allorgs)}</TableCell>
          <TableCell align="right" style={{ verticalAlign: 'top' }}>{foldOntologyData(row, "disease_info", "http://purl.obolibrary.org/obo/")}</TableCell>
          <TableCell align="right" style={{ verticalAlign: 'top' }}>{foldOntologyData(row, "ncits_info", "http://purl.obolibrary.org/obo/")}</TableCell>
          <TableCell align="right" style={{ verticalAlign: 'top' }}>{foldOntologyData(row, "go_info", "http://purl.obolibrary.org/obo/")}</TableCell>
          <TableCell align="right" style={{ verticalAlign: 'top' }}>{foldOntologyData(row, "cells_info", "http://purl.obolibrary.org/obo/")}</TableCell>
          <TableCell align="right" style={{ verticalAlign: 'top' }}>{foldOntologyData(row, "modelanats_info", "http://purl.obolibrary.org/obo/")}</TableCell>
        </TableRow>
        <TableRow>
          <TableCell style={{ paddingBottom: 0, paddingTop: 0 }} colSpan={10}>
            <Collapse in={open} timeout="auto" unmountOnExit>
              <Box sx={{ margin: 1 }}>
                <Typography variant="h6" gutterBottom component="div">
                  Details
                </Typography>

                {makeDetailsView(row)}

              </Box>
            </Collapse>
          </TableCell>
        </TableRow>
      </React.Fragment>
    );
  }
  

  function highlightedText( sentence:any, highlightAt:any)
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

      console.log(sortedHighlight);
      console.log(sentence);
      console.log(sentence.length);

      var utfChars:any = [];

      for (var i = 0; i < sentence.length; ++i)
      {
          if (sentence.charCodeAt(i) > 127)
          {
              //utfChars.push(i);
          }
      }

      console.log(utfChars)

      var utfCharsBelow = function(pos:any) {
          var ret = 0;

          for (var i = 0; i < utfChars.length; ++i)
          {
              if (utfChars[i] < pos)
              {
                  ret += 1
              }
          }

          return ret;
      }



      for (var i = 0; i < sortedHighlight.length; ++i)
      {
          var intStart = sortedHighlight[i][0];
          var intStop = sortedHighlight[i][1];
          var highColor = sortedHighlight[i][2];

          var highlightPart = sentence.substr(intStart, intStop-intStart)
          intStart = intStart - 2*utfCharsBelow(intStart);
          intStop = intStop - 2*utfCharsBelow(intStop);

          highlightPart = sentence.substr(intStart, intStop-intStart)

          var rightPart = sentence.substr(intStop, lastLeft-intStop)
          

          allParts.push(<span key={allParts.length}>{rightPart}</span>);
          allParts.push(<span style={{color:highColor, fontWeight:"bold"}} key={allParts.length}>{highlightPart}</span>);

          lastRight = intStop;
          lastLeft = intStart;
      }

      allParts.push(<span key={allParts.length}>{sentence.substr(0, sortedHighlight[sortedHighlight.length-1][0])}</span>);

      allParts.reverse();

      return <span>{allParts}</span>
  }

function prepareEvidences( allInfo:any )
  {
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
      var evStuff:any = [];

      var sortedEvs = allEvs.sort((lev:any, rev:any) => {
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

              var infoRows =prepareInfoRowPubmed(tev, evStuff.length, "");
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


              var infoRows =prepareInfoRowMirecords(tev, evStuff.length);
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

              var infoRows =prepareInfoRowMirTarBase(tev, evStuff.length);
              evStuff = evStuff.concat(infoRows);
          } else if(tev['data_source'] == "DIANA")
          {
              if (beforeEvType != 'DIANA')
              {
                  // add header
                  var headRow = <tr key={evStuff.length}>
                  <th>Gene/Organism</th>
                  <th>miRNA</th>
                  <th>Measurement Method/Type</th>
                  <th>Tissue/Cell</th>
                  <th>Link To DIANA TarBase</th>
                  </tr>;
          
                  evStuff.push(headRow);

                  beforeEvType = 'DIANA';
              }

              var infoRows = prepareInfoRowDiana(tev, evStuff.length);
              evStuff = evStuff.concat(infoRows);

          }



      }

      return <table style={{ tableLayout: "auto", width: "100%"}}>
              <tbody>
                  {evStuff}
              </tbody>
          </table>;

  }

  function prepareInfoRowMirTarBase(tev:any, idx:any)
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

  function prepareInfoRowDiana(tev:any, idx:any)
  {

      /*
              {
                  data_id: 156905
                  ​​​​​​data_source: "DIANA"
                  ​​​​​​direction: "DOWN"
                  ​​​​​​lid: "CXCR4"
                  ​​​​​​lontid: "CXCR4"
                  ​​​​​​ltype: "gene"
                  ​​​​​​measure: "DIRECT"
                  ​​​​​​method: "HITS-CLIP"
                  ​​​​​​orgs: Array [ "hsa" ]
                  ​​​​​​rid: "miR-1226-5p"
                  ​​​​​​rontid: "miR-1226-5p"
                  ​​rtype: "mirna"
                  ​​​tissue: "Kidney"
              },


                                  var headRow = <tr key={evStuff.length}>
                  <th>Gene/Organism</th>
                  <th>miRNA</th>
                  <th>Measurement Method/Type</th>
                  <th>Tissue/Cell</th>
                  </tr>;
              */

      var outRows = [];

      var longOrg = tev['orgs'].map(x => orgTLC2Long[x])

      var lowerGene = tev['lid'].toLowerCase()
      var upperGene = tev['lid'].toUpperCase()

      var dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D="+upperGene+"&genes%5B%5D="+lowerGene+"&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1"

      var infoRow = <tr key={idx+1}>
                      <td>{tev['lid']} ({tev['ltype']})<br/>{longOrg.join(", ")}</td>
                      <td>{tev['rid']} ({tev['rtype']})</td>
                      <td>{tev['method']}<br/>{tev['measure']}</td>
                      <td>
                          {
                              tev['tissue']
                          }<br/>
                          {
                              tev['cellline']
                          }
                      </td>
                      <td>
                          <a href={dianaLink} target="_blank">DIANA</a>
                      </td>
                    </tr>;

      outRows.push(infoRow);

      return outRows;
  }

function prepareInfoRowMirecords(tev:any, idx:any)
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

function prepareInfoRowPubmed(tev, idx, outlinkBase)
  {

      var outRows = [];




      var relDirection;
      if (tev['rel_direction_verb'] != null)
      {
          relDirection = tev['rel_direction_verb'];
      } else {
          relDirection = tev['rel_interaction'];
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

          orgInfo = tev['orgs'].map((d:any) => orgTLC2Long[d] || d).join(", ");
      }

      //
      var self = this;

      var infoRow = <tr key={idx}>
                      <td>{tev['rel_interaction']}, {tev['rel_category']}</td>
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
          sentRow = <tr key={idx+1}><td colSpan={4}>{highlightedText(tev['sentence'], tevHighlights)}</td></tr>;
      }

      outRows.push(sentRow);

      return outRows;

}


function makeDetailsView (row:any) {

    console.log("subcomponent");
    console.log(row);


    return (
        <div>
            <div style={{ padding: "20px" }}>

                {prepareEvidences(row)}

                </div>
        </div>
        )
}

  export default function ResultTable(input:any) {

    var data = input.data;
    console.log(data);

    return (
      <TableContainer component={Paper}>
        <Table stickyHeader aria-label="collapsible sticky table">
          <TableHead>
            <TableRow>
              <TableCell />
              <TableCell>Gene</TableCell>
              <TableCell align="right">miRNA</TableCell>
              <TableCell align="right">Document</TableCell>
              <TableCell align="right">Organisms</TableCell>
              <TableCell align="right">Disease</TableCell>
              <TableCell align="right">NCIT</TableCell>
              <TableCell align="right">Gene Ontology</TableCell>
              <TableCell align="right">Cell Ontology</TableCell>
              <TableCell align="right">Model Anatomy</TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {data.map((row:any) => (
              <Row key={row.name} row={row} />
            ))}
          </TableBody>
        </Table>
      </TableContainer>
    );
  }

/*
  const columns = useMemo(
    () => [
        {
          Header: "Found Interaction",
          columns: [
            {
              Header: "Gene/lncRNA",
              id: "lid",
              accessor: d => d.lid,
              filterMethod: (filter, rows) => 
              {
                  var elems = [ rows['lid'] ];
                  
                  var retval = matchSorter(elems, filter.value);

                  return retval.length > 0;
              },
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
              filterMethod: (filter, rows) => 
              {
                   var elems = [ rows['rid'] ];
                  
                  var retval = matchSorter(elems, filter.value);


                  return retval.length > 0;
              },
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
              },
              filterMethod: (filter, row) => {
                  var filterID = filter.id;
                  var rowData = row[filterID];

                  console.log(row);
                  console.log(filter);

                  var retval = matchSorter(rowData, filter.value);
                  console.log(retval);

                  return retval.length > 0;
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
      ],
    []
  );


  */