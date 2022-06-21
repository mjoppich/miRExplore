/* eslint no-console: "off" */

import * as path from 'path';
import { Server } from 'http';
import * as express from 'express';

import * as React from 'react';
import { renderToString } from 'react-dom/server';
import { StaticRouter as Router } from 'react-router-dom';
import { MainApp } from './pages/App';
import * as bodyparser from 'body-parser';

function toBool(input: string): boolean | undefined {
  try {
      return JSON.parse(input.toLowerCase());
  }
  catch (e) {
      return undefined;
  }
}

import * as d3 from 'd3';

// load list of genes and list of mirnas

// start the server
const port = +process.env.PORT || 65523;
const env = process.env.NODE_ENV || 'production';
const SERVER_RENDER = toBool(process.env.UNIVERSAL) || false;

console.log(port);
console.log(env);
console.log(SERVER_RENDER);

const app = express();
const server = new Server(app);

var urlEncParser = bodyparser.urlencoded({ extended: false });
var jsonParser = bodyparser.json();

// use ejs templates
app.set('view engine', 'ejs');
app.set('views', path.join(__dirname, 'views'));

// define the folder that will be used for static assets
app.use(express.static(path.join(__dirname, 'static')));

// universal routing and rendering
app.get('*', (req, res) => {
  let markup = '';
  let status = 200;

  if (SERVER_RENDER) {
    const context:any = {};
    markup = renderToString(
      <Router location={req.url} context={context}>
      <MainApp/>
      </Router>
      ,
    );
    //

    // context.url will contain the URL to redirect to if a <Redirect> was used
    if (context.url) {
      return res.redirect(302, context.url);
    }

    if (context.is404) {
      status = 404;
    }
  }

  return res.status(status).render('index', { markup });
});

server.listen(port, (err) => {
  if (err) {
    return console.error(err);
  }
  return console.info(
    `
      Server running on http://localhost:${port} [${env}]
      Universal rendering: ${SERVER_RENDER ? 'enabled' : 'disabled'}
    `);
});
