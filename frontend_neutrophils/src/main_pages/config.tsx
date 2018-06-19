

export default class config {

    static restServer = 'http://localhost'
    static restPort = '65522'

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
