/* eslint react/forbid-prop-types: "off" */

import * as React from "react"; 
import { Link } from 'react-router-dom';

export default class NotFoundPage extends React.Component<{},{}> {
  componentWillMount() {

  }

  render() {
    return (<div className="not-found">
      <h1>404</h1>
      <h2>Page not found!</h2>
      <p>
        <Link to="/">Go back to root.</Link>
      </p>
    </div>
    );
  }
};
