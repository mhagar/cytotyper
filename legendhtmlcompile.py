# This script should generate an html file out of the data fed into it
import jinja2
import pandas as pd
import sqlite3

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "template.html.jinja"
template = templateEnv.get_template(TEMPLATE_FILE)

colors = ['ff0000', 'e3256b', 'c4717a', '09cb95', '65b0bf', 'b87333', 'cf6ae6']

conn = sqlite3.connect('StzF-SSN-Nov-26-2020.sqlite')

sqlquery = "family, seq_len"


def average_size(pfam):
    size_df = pd.read_sql_query(f"SELECT {sqlquery} FROM neighbors \
                                WHERE family LIKE '%{pfam}%'\
                                UNION ALL SELECT {sqlquery} FROM attributes \
                                WHERE family LIKE '%{pfam}%' ",
                                conn)
    return int(size_df['seq_len'].mean())


class Enzyme:
    def __init__(self, label, pfam, color):
        self.label = label
        self.pfam = pfam
        self.color = colors[color]
        self.size = average_size(pfam[0])


def generate_enzyme_colors(queryfams, queryfam_labels):
    output = {}
    for key, label in queryfam_labels.items():
        output[key] = Enzyme(label, queryfams[key], key)

    # import pdb; pdb.set_trace()
    return output


def generate_bgc_legend(genotypes, genotype_count):
    output = {}
    for key, item in genotype_count.items():
        output[key] = [genotypes[key], genotype_count[key]]
    return output


def export_data_report(filename, queryfams, queryfam_labels,
                       genotypes, genotype_count):
    templateinputs = {
     'enzymes': generate_enzyme_colors(queryfams, queryfam_labels),
     'genotypes': generate_bgc_legend(genotypes, genotype_count),
    }
    outputText = template.render(templateinputs)
    f = open(filename, 'w')
    final = outputText
    f.write(final)
    f.close
