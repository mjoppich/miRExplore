import os

from flask import Flask, request
# set the project root directory as the static folder, you can set others.

dataurl = '/home/mjoppich/git/miRExplore/frontend_neutrophils/src/static/'

app = Flask(__name__, static_folder=dataurl, static_url_path='/static')

@app.route('/')
def root():

    retFile = 'index.html'

    return app.send_static_file(retFile)

app.run()