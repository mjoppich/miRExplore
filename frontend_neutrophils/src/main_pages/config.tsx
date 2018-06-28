

export default class config {

    static restServer = 'http://localhost'
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
