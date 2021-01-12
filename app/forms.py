from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired
from wtforms import BooleanField, SubmitField
from wtforms.validators import DataRequired


class TestForm(FlaskForm):
    user_upload = FileField(validators=[FileRequired()])
    match_direction = BooleanField('Match direction?')
    submit = SubmitField('GO!!')
