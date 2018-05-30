

export default class config {

    static restServer = 'http://bioclient1.bio.ifi.lmu.de'
    static restPort = '65526'

    static getRestAddress()
    {
        return this.restServer + ":" + this.restPort ;
    }

    static axiosConfig = {
        crossdomain: true,
        headers: {
            "Access-Control-Allow-Origin": "*"
        }
      };
  
}
