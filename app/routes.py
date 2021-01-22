from flask import render_template, redirect, send_file
from app import app
from app.forms import TestForm
import os.path
import filewrap
import cytotyper_test
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
    rawdata = form.user_upload.data
    if form.validate_on_submit() and allowed_upload(rawdata.filename):
        # import pdb; pdb.set_trace()
        cytotyper_test.header(rawdata).to_csv('app/output/output.csv')
        # rawdata.save(os.path.join(app.config['UPLOAD_DIR'], rawdata.filename))
        filewrap.rmfile('app/output.zip')
        filewrap.zip_wrap('app/output')
        return send_file('output.zip', attachment_filename='results.zip')
    return render_template('index.html',
                           title='Home',
                           form=form)
