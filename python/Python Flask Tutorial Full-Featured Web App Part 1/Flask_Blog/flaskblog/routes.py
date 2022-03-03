from flask import render_template, url_for, flash, redirect
from flaskblog.forms import RegistrationForm, LoginForm
from flaskblog.forms import User
from flaskblog import app, db, bcrypt


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
    return render_template('about.html', title='About')


@app.route('/register', methods=['GET', 'POST'])
def register():
    form = RegistrationForm()
    if form.validate_on_submit():
        hashed_password = bcrypt.generate_password_hash(form.password.data).decode('utf-8')
        user = User(username=form.username.data, email=form.email.data, password=hashed_password)
        db.session.add(user)
        db.session.commit()
        (flash(f'Your account has been created. You are now able to login!', 'success'))
        return redirect(url_for('login'))
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
