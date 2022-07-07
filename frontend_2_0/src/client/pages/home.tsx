import * as React from 'react'

import Container from '@mui/material/Container'
import Typography from '@mui/material/Typography'
import Box from '@mui/material/Box'

import QueryComponent from '../components/QueryComponent'

export interface HomeProps { };
export interface HomeState { allQueries: any};

class Home extends React.Component<HomeProps, HomeState> {

  constructor(props: any) {
    super(props);

    var self=this;
    this.state = {allQueries: [<QueryComponent key={0}/>]};
  }


  newQuery()
  {
    var self=this;

      this.state.allQueries.push(<QueryComponent key={self.state.allQueries.length}/>);
  }

  clearQueries()
  {
      this.setState({allQueries: []});
  }


  render() {

    return (<div style={{padding:40}}>
    <Container maxWidth="false">
      <Box sx={{ my: 4 }}>
        <Typography variant="h2" component="h2" gutterBottom>
        miRExplore 🧪🧬
        </Typography>

        <Typography variant="body1" gutterBottom>
        <div>
                <span>
                    <h4>General Information</h4>
                    <p>
                        miRExplore is an interactive web tool, which allows the user to browse miRNA - gene interactions in human and mouse. 
                    </p>
                </span>
                <span>
                    <h4>Included Databases</h4>
                    <ul>
                        <li><b>miRBase</b>: 1918 (human), 1227 (mouse)</li>
                        <li><b>miRecords</b></li>
                        <li><b>PubMed and PubMed Central</b></li>
                    </ul>
                </span>
          </div>
          </Typography>

      </Box>
    </Container>

    <Container maxWidth="false">
      <Box sx={{ my: 4 }}>
        <Typography variant="h2" component="h2" gutterBottom>
          Analysis Builder
        </Typography>

        <Typography variant="body1" gutterBottom>
        <p>Using a query you can query our database for interactions.</p>
        <p>You must select an entity for which you would like to explore interactions. This might be one term, but can also be multiple terms.</p>
        <p>You can restrict your search to such evidences, which must be for a specific organism, or must include a specific messenger or category.</p>
        <p>You can show a interaction network for your selected evidences, or a sankey chart, showing you also how the interactions are connected to messengers or categories</p>
        <p>Enjoy exploring your miRNAs!</p>
        </Typography>
      </Box>
    </Container>

    <Container maxWidth="false">
      <Box sx={{ my: 4 }}>
        <Typography variant="h2" component="h2" gutterBottom>
          Queries
        </Typography>

        <Typography variant="body1" gutterBottom>
        {this.state.allQueries}
        </Typography>
      </Box>
    </Container>
    </div>
  )

  }

}



export default Home
