const path = require('path')
const webpack = require('webpack');
const HtmlWebpackPlugin = require('html-webpack-plugin')
const MiniCssExtractPlugin = require('mini-css-extract-plugin')

const outputDirectory = 'dist'

module.exports = {
  entry: ['./src/client/index.tsx'],
  output: {
    path: path.join(__dirname, outputDirectory),
    filename: './js/[name].bundle.js',
  },
  devtool: 'source-map',
  module: {
    rules: [
      {
        test: /\.(js|jsx|ts|tsx)$/,
        exclude: /node_modules/,
        use: {
          loader: 'babel-loader',
        },
      },
      {
        enforce: 'pre',
        test: /\.js$/,
        loader: 'source-map-loader',
      },
      {
        test: /\.css$/,
        use: [{ loader: 'style-loader' }, { loader: 'css-loader' }],
      },
    ],
  },
  resolve: {
    extensions: ['*', '.ts', '.tsx', '.js', '.jsx', '.json'],
  },
  watchOptions: {
    aggregateTimeout: 300,
    poll: 2000,
    ignored: /node_modules/,
  },
  devServer: {
    port: 8000,
    open: true,
    hot: true,
    proxy: {
      '/api': {
        target: 'http://localhost:8001',
        secure: false,
        changeOrigin: true,
      },
    },
  },
  plugins: [
    new webpack.DefinePlugin({
      PUBLIC_URL: '/',
      NODE_ENV: process.env.NODE_ENV,
    }),
    new HtmlWebpackPlugin({
      inject: true,
      template: './public/index.html',
      favicon: './public/favicons/favicon.ico',
      title: 'express-typescript-react',
    }),
    new MiniCssExtractPlugin({
      filename: './css/[name].css',
      chunkFilename: './css/[id].css',
    }),
  ],
}
