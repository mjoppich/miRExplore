

export default class config {

    static restServer = 'http://localhost'//'https://turingwww.bio.ifi.lmu.de'
    static restPort = '65522'
    static restFolder = '';//'neutrophils'

    static getRestAddress()
    {
        var basePart = this.restServer;

        if ((this.restPort != null) && (this.restPort.length > 0))
        {
                basePart += ":" + this.restPort;
        }

        if ((this.restFolder != null) && (this.restFolder.length > 0))
        {
            basePart += "/" + this.restFolder;
        }

        return basePart;
    }

    static axiosConfig = {
        crossdomain: true,
        headers: {
            "Access-Control-Allow-Origin": "*"
        }
      };

}