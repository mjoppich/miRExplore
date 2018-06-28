

export default class config {

    static restServer = 'http://localhost'//'http://bioclient1.bio.ifi.lmu.de'
    static restPort = '65522'
    static restFolder = null;//'neutrophils'

    static getRestAddress()
    {
        var basePart = this.restServer + ":" + this.restPort;

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
