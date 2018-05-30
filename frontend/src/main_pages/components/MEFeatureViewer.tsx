import * as React from "react"; 
import * as ReactDOM from "react-dom";
import FeatureViewer from './FeatureViewer';


export interface MEFeatureViewerProps { features: any, location:any}
export interface MEFeatureViewerState { }



export default class MEFeatureViewer extends React.Component<MEFeatureViewerProps, MEFeatureViewerState> {

    constructor(props) {
        super(props);

        console.log(this.props);

    }

    randomSeq(seqLen: number) {
        var text = "";
        var possible = "ABCDEFGHIJKLMNOPQRSTUVWXY";
      
        for (var i = 0; i < seqLen; i++)
          text += possible.charAt(Math.floor(Math.random() * possible.length));
      
        return text;
      }
      

    componentDidMount()
    {
        const de = this.refs.fv1;

        console.log(de);
        console.log($("#fv1"));

        var ft = new FeatureViewer(this.randomSeq(200), "#fv1", {
             showAxis: true,
             showSequence: true,
             brushActive: true, //zoom
             toolbar:true, //current zoom & mouse position
             bubbleHelp:true, 
             zoomMax:50 //define the maximum range of the zoom
         });

         ft.addFeature({
            data: [{x:20,y:32},{x:46,y:100},{x:123,y:167}],
            name: "let-7 binding sites",
            className: "test1",
            color: "#0F8292",
            type: "rect",
            filter: "type1"
        });

        ft.addFeature({
            data: [{x:15,y:23},{x:67,y:75},{x:99,y:130}],
            name: "miR-130 binding sites",
            className: "test1",
            color: "#378a72",
            type: "rect",
            filter: "type1"
        });
        

    }

    render() {

        console.log("ME FeatureViewer Render");

        return (<div><div id="fv1" style={{width: "1000px"}}></div></div>);
    }

}