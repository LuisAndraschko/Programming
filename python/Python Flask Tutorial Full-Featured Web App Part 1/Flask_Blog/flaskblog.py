from turtle import title
from flask import Flask, render_template, url_for, flash, redirect
from forms import RegistrationForm, LoginForm
import os

app = Flask(__name__)

SECRET_KEY = os.urandom(32)
app.config['SECRET_KEY'] = SECRET_KEY

posts = [
        {
            'author': 'Luís Andraschko',
            'title': 'Blog Post 1',
            'content': 'First post content',
            'date_posted': 'April 20, 2018'
        },
        {
            'author': 'Corey Schafer',
            'title': 'Blog Post 2',
            'content': 'Second post content',
            'date_posted': 'April 21, 2018'
        }
    ]


@app.route("/")
@app.route("/home")
def home():
    return render_template('home.html', posts=posts)

@app.route("/about")
def about():
    return render_template('about.html', 
                           title='About')

@app.route('/register', methods=['GET', 'POST'])
def register():
    form = RegistrationForm()
    if form.validate_on_submit():
        (flash(f'Account created for {form.username.data}!', 'success'))
        return redirect(url_for('home')) 
        
    return render_template('register.html', 
                           title='Register', form=form)

@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        if form.email.data == 'admin@gmail.com' and form.password.data == 'password':
            (flash(f'Loggeg in as {form.email.data}.', 'success'))
            return redirect(url_for('home'))
        else:
            (flash(f'Could not find a match, try again.', 'danger'))
    return render_template('login.html', 
                           title='Login', form=form)


if __name__ == '__main__':
    app.run(debug=True) 