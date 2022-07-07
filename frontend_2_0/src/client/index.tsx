import 'core-js/stable'
import 'regenerator-runtime/runtime'

import * as React from 'react'
import * as ReactDOMClient from 'react-dom/client'
import './styles/stylesheet.css'

import App from './containers/app'

const container = document.getElementById('root') as Element


// Create a root.
const root = ReactDOMClient.createRoot(container)


root.render(<App />)
