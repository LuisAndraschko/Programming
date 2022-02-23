"""
CMD commands to start application:
[1]
cd /D D:\Estudos\Programming\python\Python Flask Tutorial Full-Featured Web App Part 1\Flask_Blog\ 

[2]
set FLASK_APP=flaskblog.py

[3]
flask run
"""
"""
To open application on browser use the URL:

http://localhost:5000/
"""
"""
To refresh a browser page it's needed to reload the webserver via command line:

CRTL+C

redo start application commands
"""
"""
In order to avoid all this trouble whenever I want to reload the browser page we need to set another environment variable called FLASK_DEBUG equal to 1 instead of FLASK_APP="py file name"[2] and run flask[3] we can simply reload the browser page whenever we want to see a change
"""
"""
In order to accelerate this process we can set the following if statement:

if __name__ == '__main__':
    app.run(debug=True) 

in this file and whenever we want to start the application we simply type on the CMD:

python flaskblog.py
"""
from flask import Flask

app = Flask(__name__)

@app.route("/")
@app.route("/home")
def homepage():
    return "<h1>This is your home page</h1>"

@app.route("/about")
def about():
    return "<p>This is your about page</p>"

if __name__ == '__main__':
    app.run(debug=True) 