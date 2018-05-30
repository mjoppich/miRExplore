/* global window document */

import * as React from 'react';
import { render } from 'react-dom';
import { BrowserRouter as Router } from 'react-router-dom';
import { MainApp } from './pages/App';

import $ from 'jquery'
import boostrap from 'bootstrap'


const AppClient = () => 
<Router>
  <MainApp/>
</Router>
;

window.onload = () => {
  render(<AppClient />, document.getElementById('main'));
};
