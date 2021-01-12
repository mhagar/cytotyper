from flask import render_template, redirect
from app import app
from app.forms import TestForm
import os.path
from werkzeug.utils import secure_filename


def allowed_upload(filename):
    if '.' in filename:
        if filename.rsplit('.', 1)[1].lower() == 'sqlite':
            return True
    else:
        return False


@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = TestForm()
    f = form.user_upload.data
    if form.validate_on_submit() and allowed_upload(f.filename):
        # import pdb; pdb.set_trace()
        f.save(os.path.join(app.config['UPLOAD_DIR'], f.filename))
        return redirect('/index')
    return render_template('index.html',
                           title='Home',
                           form=form)
