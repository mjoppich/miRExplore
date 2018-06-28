

export default class config {

    static restServer = 'http://bioclient1.bio.ifi.lmu.de'
    static restPort = '80'
    static restFolder = 'neutrophils'

    static getRestAddress()
    {
        return this.restServer + ":" + this.restPort + "/" + this.restFolder + "/";
    }

    static axiosConfig = {
        crossdomain: true,
        headers: {
            "Access-Control-Allow-Origin": "*"
        }
      };
  
}
