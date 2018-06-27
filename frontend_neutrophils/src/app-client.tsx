/* global window document */

import * as React from 'react';
import { render } from 'react-dom';
import { BrowserRouter as Router } from 'react-router-dom';
import { MainApp } from './pages/App';

import "jquery";
import "bootstrap";
import 'bootstrap/dist/css/bootstrap.min.css';


const AppClient = () => 
<Router basename={'/neutrophils'}>
  <MainApp/>
</Router>
;

window.onload = () => {
  render(<AppClient />, document.getElementById('main'));
};
