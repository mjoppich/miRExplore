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
const port = +process.env.PORT || 3000;
const env = process.env.NODE_ENV || 'production';
const SERVER_RENDER = toBool(process.env.UNIVERSAL) || false;

console.log(port);
console.log(env);
console.log(SERVER_RENDER);

const app = express();
const server = new Server(app);

var urlEncParser = bodyparser.urlencoded({ extended: false });
var jsonParser = bodyparser.json();

import * as n4j from 'neo4j-driver';
var n4jdriver = n4j.v1.driver("bolt://localhost", n4j.v1.auth.basic("neo4j", "neo4j2"));
var _ = require('lodash');

// use ejs templates
app.set('view engine', 'ejs');
app.set('views', path.join(__dirname, 'views'));

// define the folder that will be used for static assets
app.use(express.static(path.join(__dirname, 'static')));


app.get("/api", (req, res) => {
  const session = n4jdriver.session();

  session.run("MATCH (n) RETURN distinct labels(n) as desc, count(*) as num;").then(
    result => {
      session.close();

      var retValues = {};
    
      for (let elem of result.records)
      {
        console.log(elem);

        for (let descr of elem.get('desc'))
        {

          if (!(descr in retValues))
          {
            retValues[descr] = 0
          }

          retValues[ descr ] += n4j.v1.integer.toNumber(elem.get('num'));

        }
      }
    
      // on application exit:
      n4jdriver.close();

      res.send(retValues);

  });  

});

app.all('/api/acinput', urlEncParser, (req, res, next) => {

  console.log(req.params);
  console.log(req.body);
  console.log(req.query);

  res.send({
    params: req.params,
    bodyparams: req.body,
    query: req.query
  });
  
});

app.get("/api/query", (req, res) => {
  const session = n4jdriver.session();

  session.run(
    'MATCH (n:GENE)-[]-(m:EVIDENCE) RETURN n, collect(m) as evidences LIMIT {limit}', {limit: 10})
    .then(results => {
      session.close();

      var nodes = new Array<{id: any, title:string, label:string, r: number, group: number}>();
      var rels = new Array();

      //res.send(results);

      results.records.forEach(res => {
        nodes.push({id: res.get('n').properties.id, title: res.get('n').name, label: 'gene', 'r': 8, 'group': 1});
        var target = nodes.length-1;

        res.get('evidences').forEach(name => {
          console.log(name);

          var sourceIdx = [];

          var actor = {id: name.properties.id, title: name, label: 'evidence', 'r': 8, 'group': 1};
          var sourceValues = nodes.filter((elem, idx) => {

            if (elem.id == name.properties.id)
            {
              sourceIdx.push(idx);
            }
        });
          var source = -1;

          if (sourceIdx.length == 0) {
            nodes.push(actor);
            source = nodes.length-1;

            rels.push({source: source, target: target, d: 10})

          } else {
            for (let idx of sourceIdx)
            {
              rels.push({source: idx, target: target, d: 10})
            }
          }
        })
      });

      console.log("Nodes: " + nodes.length);

      var result = {nodes: nodes, links: rels};
      res.send(result);
});

});

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
