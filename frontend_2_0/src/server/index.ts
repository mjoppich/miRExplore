import express, { json, Request, Response, Router, Express } from 'express'

// setting up express server
const app: Express = express()

app.use(json())

const port: number = Number(process.env.PORT) || 8001;
const PROD = process.env.NODE_ENV === 'production'


if (PROD) {
  app.use('/', express.static('dist'))
}

app.get('/api/message', (req, res) => {
  res.send('Build something amazing! 🚀')
})

app.listen(port, () => {
    console.log(`\n\n🚀  Server running at http://localhost:${port}.\n\n`)
})