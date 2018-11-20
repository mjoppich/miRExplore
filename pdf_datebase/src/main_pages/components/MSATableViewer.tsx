import * as React from "react"; 
import FlatButton  from 'material-ui/FlatButton';
import DownloadArchive from 'material-ui/svg-icons/content/archive';
import Toggle from 'material-ui/Toggle';

export interface DownloadButtonProps { filename: string, getDownloadContent: any }


class DownloadButton extends React.Component<DownloadButtonProps, {}> {

    makeAlignmentDownload()
    {
        var fileContents = this.props.getDownloadContent();

        var filetype = "text/plain";

        var a = document.createElement("a");
        var dataURI = "data:" + filetype +
            ";base64," + btoa(fileContents);
        a.href = dataURI;
        a['download'] = this.props.filename;
        var e = document.createEvent("MouseEvents");
        // Use of deprecated function to satisfy TypeScript.
        e.initMouseEvent("click", true, false,
            document.defaultView, 0, 0, 0, 0, 0,
            false, false, false, false, 0, null);
        a.dispatchEvent(e);

        a.parentNode.removeChild(a);
    }

    render()
    {
        return <FlatButton 
        label=""
        labelPosition="before"
        primary={true}
        icon={<DownloadArchive />}
        onClick={() => this.makeAlignmentDownload()}
        />
        
    }

}

export interface MSATAbleViewerProps { alignments: any }
export interface MSATAbleViewerState { showAAAlign: boolean, showNTAlign:boolean }


export default class MSATableViewer extends React.Component<MSATAbleViewerProps, MSATAbleViewerState> {

    xrefcaptions: any;
    goNS2Name: any;

    constructor(props: MSATAbleViewerProps) {
        super(props);

        this.xrefcaptions = {
            'GO': 'Gene Ontology (GO)',
            'TMGO': 'Textmining GO',

            'Pfam': 'Protein Families',
            'Interpro': 'Interpro'
        };

        this.goNS2Name = {
            'biological_process': 'Biological Process',
            'molecular_function': 'Molecular Funciton',
            'cellular_component': 'Cellular Component'
        }
    }

    getGONS2Name( namespaceID, defaultVal )
    {
        if (this.goNS2Name[namespaceID] !== undefined)
        {
            return this.goNS2Name[namespaceID];
        }

        return defaultVal;
        
    }


    componentDidMount()
    {
    }

    componentWillMount()
    {
        this.setState({showAAAlign: true, showNTAlign: false});
    }

    makeGOLinks( rowAlign, idCat )
    {

        if (!(idCat in rowAlign.xrefs))
        {
            return <div></div>;
        }

        var allGOClasses = Object.keys(rowAlign.xrefs[idCat])

        console.log("all GO Classes");
        console.log(allGOClasses);
        console.log(allGOClasses.map((namespaceID, nsi) => rowAlign.xrefs[idCat][namespaceID]));

        var self = this;


        var elemsByNamespace = allGOClasses.map((namespaceID, nsi) => 
            <div key={nsi}>
                <h3 style={{margin: 0}}>{self.getGONS2Name(namespaceID, namespaceID)}</h3>
                {rowAlign.xrefs[idCat][namespaceID].map(
                    (x, i) => 
                                <div key={i}>
                                    <a href={"http://amigo.geneontology.org/amigo/term/"+x}>{x}</a><br/>
                                    <span>{rowAlign.xrefs["GOTERMS"][x]}</span><br/>
                                </div>
                )}
            </div>
        );

        return <div>{elemsByNamespace}</div>
    }

    printMSA()
    {
        var outStr = "";
        for (var i = 0; i < this.props.alignments.msa.length; ++i)
        {
            var rowAlign = this.props.alignments.msa[i];

            var rowID = rowAlign.entryID;
            var alignSeq = rowAlign.alignment;

            outStr += ">" + rowID + "\n" + alignSeq + "\n";
        }

        return outStr;
    }

    printSeq( seqid )
    {
        var outStr = "";
        for (var i = 0; i < this.props.alignments.msa.length; ++i)
        {
            var rowAlign = this.props.alignments.msa[i];

            var rowID = rowAlign.entryID;
            var alignSeq = rowAlign.alignment;

            if (rowID == seqid)
            {
                outStr += ">" + rowID + "\n" + alignSeq + "\n";
            }
        }

        return outStr;
    }

    makeProtIDInfo(rowAlign)
    {

        let thisRowID = rowAlign.entryID;

        var recID = rowAlign.recordID;
        var recStart = rowAlign.genomicStart;
        var recEnd = rowAlign.genomicEnd;
        var recStrand = "?";
        
        if (rowAlign.genomicStrand == -1)
        {
            recStrand = "-";
        } else if (rowAlign.genomicStrand == 1)
        {
            recStrand = "+"
        }

        var posInfo = recID + " " + recStart + "-" + recEnd + ":" + recStrand;

        var names = [];

        for (var i = 0; i < rowAlign.entryNames.length; ++i)
        {
            names.push(<li key={i}>rowAlign.entryNames[i]</li>);
        }

        return <div>
                    <p style={{margin: 0}}>{rowAlign.entryID}</p>
                    <DownloadButton filename={thisRowID + ".fa"} getDownloadContent={() => this.printSeq(thisRowID)}/>
                    <p>{posInfo}</p>
                    <p>
                        <ul>{names}</ul>
                    </p>
                </div>
    }

    makeTSSTable(tss)
    {
        tss = tss.sort(function(a, b) {
            return a.TSSID > b.TSSID;
        })

        /*

            {
                "LOCUSTAG": "HP_1523",
                "PROPS": [
                    "secondary",
                    "enriched"
                ],
                "SEQ": "GGCTTTTTCATCTTCTTCTTCATGCTCAATTTTTCTTATATCATTATTCGC",
                "STRAND": "-",
                "TSS": 1604179,
                "TSSID": "TSS2401"
            }
        */

        var tssStuff = [];
        for (var i = 0; i < tss.length; ++i)
        {
            var tssElem = tss[i];

            tssStuff.push(
                <tr key={i}>
                    <td style={{verticalAlign: 'top'}}>{tssElem.TSSID}</td>
                    <td style={{verticalAlign: 'top'}}>{tssElem.LOCUSTAG}</td>
                    <td style={{verticalAlign: 'top'}}>{tssElem.TSS}</td>
                    <td style={{verticalAlign: 'top'}}>{tssElem.STRAND}</td>
                    <td style={{verticalAlign: 'top'}}><pre style={{margin: "0"}}>{tssElem.SEQ}</pre></td>
                    <td style={{verticalAlign: 'top'}}>{tssElem.PROPS.join(', ')}</td>

                </tr>
            );
        }

        return <table style={{ tableLayout: "auto", width: "100%"}}>
                <tbody>
                    <tr style={{textAlign: 'left'}}>
                        <th>TSS-ID</th>
                        <th>Locus-Tag</th>
                        <th>TSS</th>
                        <th>Strand</th>
                        <th>Sequence -50 nt upstream + TSS (51nt)</th>
                        <th>Properties</th>
                    </tr>
                    {tssStuff}
                </tbody>
            </table>;
    }

    makeOperonsTable(operons)
    {
        var opInfo = operons;

        // sort operons
        opInfo = opInfo.sort(function(a, b) {
            return a.OPERONID > b.OPERONID;
        })

        /*

        {
                "DOOR": [],
                "GENES": [
                    "HP_1519",
                    "HP_1520",
                    "HP_1521",
                    "HP_1522",
                    "HP_1523"
                ],
                "OPERONID": "OPERON505",
                "comment": "",
                "evaluation": "new",
                "strand": "-"
        }

        */

        var opStuff = [];
        for (var i = 0; i < opInfo.length; ++i)
        {
            var opElem = opInfo[i];

            opStuff.push(
                <tr key={i}>
                    <td style={{verticalAlign: 'top'}}>{opElem.OPERONID}</td>
                    <td style={{verticalAlign: 'top'}}>{opElem.GENES.map((x, i) => x).join(", ")}</td>
                    <td style={{verticalAlign: 'top'}}>{opElem.strand}</td>
                    <td style={{verticalAlign: 'top'}}>{opElem.evaluation}</td>
                    <td style={{verticalAlign: 'top'}}>{opElem.comment}</td>
                    <td style={{verticalAlign: 'top'}}>{opElem.DOOR}</td>
                </tr>
            );
        }

        return <table style={{ tableLayout: "fixed", width: "100%"}}>
                <tbody>
                    <tr style={{textAlign: 'left'}}>
                        <th>OperonID</th>
                        <th>Genes</th>
                        <th>Strand</th>
                        <th>Evaluation</th>
                        <th>Comment</th>
                        <th>DOOR ID</th>
                    </tr>
                    {opStuff}
                </tbody>
            </table>;
    }

    render() {

        var TSSinfo = <p>No TSS Info available</p>;
        var operonInfo = <p>No Operon Info available</p>;

        if (this.props.alignments.TSS)
        {
            var tssStuff = this.makeTSSTable(this.props.alignments.TSS);

            TSSinfo = <div>{tssStuff}</div>;
        }

        if (this.props.alignments.OPERONS)
        {
            var opStuff = this.makeOperonsTable(this.props.alignments.OPERONS);

            operonInfo = <div>{opStuff}</div>;
        }



        var allRows = []
        var allXrefs = []

        allXrefs.push(<tr key={0} style={{textAlign: 'left'}}>
        <th style={{width: "100px"}}>Organism</th>
        <th style={{width: "150px"}}>Gene/Protein Name</th>
        <th style={{width: "100px"}}>UniProt</th>

        <th>{this.xrefcaptions.GO}</th>
        <th>{this.xrefcaptions.TMGO}</th>

        <th>{this.xrefcaptions.Pfam}</th>
        <th>{this.xrefcaptions.Interpro}</th>
        </tr>)

        for (var i = 0; i < this.props.alignments.msa.length; ++i)
        {
            var rowAlign = this.props.alignments.msa[i];

            var aaSEQPart = <pre style={{display: "none"}}></pre>;
            var ntSEQPart = <pre style={{display: "none"}}></pre>;

            if (this.state.showAAAlign)
            {
                aaSEQPart = <pre style={{margin: "0"}}>{rowAlign.alignment}</pre>;
            }
            if (this.state.showNTAlign)
            {
                ntSEQPart = <pre style={{margin: "0"}}>{rowAlign.alignmentNT[0]}</pre>;

                if (this.state.showAAAlign)
                {
                    aaSEQPart = <pre style={{margin: "0"}}>{rowAlign.alignmentNT[1]}</pre>;
                }

            }
            
            allRows.push(
            <tr key={i}>
                <th style={{ position:"absolute", left:0, width: "100px", verticalAlign: "top", borderTop: "1px solid #ccc"}}>{rowAlign.entryID}</th>
                <td style={{ background:"#e7e7d7", verticalAlign: "top", borderTop: "1px solid #ccc"}}>
                {aaSEQPart}
                {ntSEQPart}
                </td>
            </tr>
            );

            var xKeys = Object.keys(rowAlign.xrefs);

            var xrefvals = [];

            for (var l=0; l < xKeys.length; ++l)
            {
                var xkey = xKeys[l];
                var xKeyVals = rowAlign.xrefs[xkey];

                for (var xv=0; xv < rowAlign.xrefs[xkey].length; ++xv)
                {
                    xrefvals.push( rowAlign.xrefs[xkey][xv] )
                }
            }


            allXrefs.push(<tr key={allXrefs.length}>
                <td style={{verticalAlign: 'top'}} ><p style={{margin: 0}}>{rowAlign.organismName}</p><p>{rowAlign.organismID}</p></td>
                <td style={{verticalAlign: 'top'}}>{this.makeProtIDInfo(rowAlign)}</td>
                <td style={{verticalAlign: 'top'}}>{rowAlign.xrefs['Uniprot'].map((x, i) => <div key={i}><a href={"http://www.uniprot.org/uniprot/"+x}>{x}</a><br/></div>)}</td>

                <td style={{verticalAlign: 'top'}}>{this.makeGOLinks(rowAlign, 'GO')}</td>
                <td style={{verticalAlign: 'top'}}>{this.makeGOLinks(rowAlign, 'TMGO')}</td>

                <td style={{verticalAlign: 'top'}}>{rowAlign.xrefs['Pfam'].map((x, i) => <div key={i}><a href={"https://pfam.xfam.org/family/"+x}>{x}</a><br/></div>)}</td>
                <td style={{verticalAlign: 'top'}}>{rowAlign.xrefs['Interpro'].map((x, i) => <div key={i}><a href={"https://www.ebi.ac.uk/interpro/entry/"+x}>{x}</a><br/></div>)}</td>

            </tr>);

        }

        return (
            <div>
            <h1>{this.props.alignments.homid}
            <DownloadButton filename={this.props.alignments.homid + ".fa"} getDownloadContent={() => this.printMSA()}/>
                </h1>
            <div>
            <Toggle
                label="Show AA Sequence"
                defaultToggled={true}
                toggled={this.state.showAAAlign}
                onToggle={(event, newValue) => this.setState({showAAAlign: newValue || (!newValue && !this.state.showNTAlign)})}
                />
            <Toggle
                label="Show NT Sequence"
                defaultToggled={false}
                toggled={this.state.showNTAlign}
                onToggle={(event, newValue) => this.setState({showAAAlign: this.state.showAAAlign || (!this.state.showAAAlign && !newValue), showNTAlign: newValue})}

/>
            </div>
            <div style={{position:"relative"}}>
                <div style={{overflowX:"scroll", overflowY:"visible", width:"90%", marginLeft:"110px"}}>
                    <table style={{ tableLayout: "auto", width: "100%"}}><tbody>
                        {allRows}
                        </tbody>
                    </table>
                </div>
        </div>
        <div>
            <h3>Cross-References</h3>
        <table style={{ tableLayout: "fixed", width: "100%"}}><tbody>
                        {allXrefs}
                        </tbody>
                    </table>

            <h3>Operon Information</h3>
            {operonInfo}
            <h3>Transcription Start Sites</h3>
            {TSSinfo}
            </div>
        </div>);
    }
}