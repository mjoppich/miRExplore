import * as React from 'react'

import Container from '@mui/material/Container'
import Typography from '@mui/material/Typography'
import Box from '@mui/material/Box'


class HelpPage extends React.Component<{ }, {}> {

  constructor(props: any) {
    super(props);

    var self=this;
  }


  newQuery()
  {
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
          miRExplore
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

    </div>
  )

  }

}



export default HelpPage
