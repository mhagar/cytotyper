import pandas as pd
import sqlite3
import legendhtmlcompile

# for debugging purposes in terminal
pd.set_option('display.max_columns', 500)

# connects to the SQLite file
conn = sqlite3.connect('StzF-SSN-Nov-26-2020.sqlite')

# generates a pandas dataframe from the 'neighbors' table in the sqlite DB, and
# extracting the following columns from that table:
sqlquery = "id, num, family, ipro_family, rel_start, direction, gene_key"
test_df = pd.read_sql_query(f"SELECT {sqlquery} FROM neighbors", conn)

# so it turns out, in the 'attributes' table in the sqlite db, the index
# corresponds to the gene_key! (i.e. gene_key 573, which refers to the
# kutzneria BGC, the StzF-homolog is found in table 'attributes' row 573.)
sqlquery = "sort_key, accession, num, family, ipro_family, rel_start, direction"
organism_df = pd.read_sql_query(f"SELECT {sqlquery} FROM attributes", conn,
                                index_col='sort_key')


# this function takes 'AHMAD-BORSCHT-CAR' and returns [AHMAD, BORSCHT, CAR]
def split_dash(x):
    x = x.split('-')
    return x


# this line uses the above function to convert the 'family' and 'ipro_family'
# columns in the dataframe into list objects
test_df['family'] = test_df['family'].apply(split_dash)
test_df['ipro_family'] = test_df['ipro_family'].apply(split_dash)

# this line takes everything in "ipro_family" and also adds it to "family"
# then gets rid of the "ipro_family" column
# (makes it easier to search later)
test_df['family'] += test_df['ipro_family']
test_df.drop(labels='ipro_family', axis=1, inplace=True)


# this next block opens 'query.csv' and retrieves a table formatted like this:
#   Label   |   Pfam codes (separated by a ',')
# You can input multiple Pfam codes as 'synonyms' that you want to be grouped
# under one label.
querycsv = pd.read_csv('query.csv')
queryfam_labels = {}
queryfams = []

for index, item in enumerate(querycsv['Label'].tolist()):
    queryfam_labels[index] = item
    pass

for item in querycsv['Pfams'].tolist():
    queryfams.append(item.split(','))
    pass
# at the end of this block, you're left with QUERYFAMS_LABELS as a dict where
# the format is: { 0: Label, etc..}
# QUERYFAMS is a list of lists, such that synonyms are grouped together like:
# [[PfamA, PfamB, PfamC], [PfamD, PfamE], [PfamF]]
# The index corresponds to the 'label' in QUERYFAMS_LABELS


# this next block is the meat and potatoes.

gene_key_hits = {}
genotypes = []
genotype_count = {}
genotype_genes = {}
cytotable_genotypes = {}

# go thru each row in the df..
# for testing: kutz is (18724, 18763)
for i in range(0, len(test_df)):
    # hit = False
    # for every set of synonyms in the query..
    for queryfamkey, queryfamsynonyms in enumerate(queryfams):
        # this line checks if anything in a given synonym set matches the
        # 'family' cell in that row..
        hits = [x for x in queryfamsynonyms if x in test_df.loc[i, 'family']]
        if hits:
            direction = test_df.loc[i, 'direction']
            gene_key = test_df.loc[i, 'gene_key']
            rel_start = test_df.loc[i, 'rel_start']
            # if the direction of a gene is the same as StzF's..
            if (direction == organism_df.loc[gene_key, 'direction'] and
               abs((test_df.loc[i, 'num'] - organism_df.loc[gene_key, 'num'] < 10))):
                if gene_key not in gene_key_hits:
                    gene_key_hits[gene_key] = []
                    # (this next bit is for later..)
                    cytotable_genotypes[gene_key] = [organism_df.loc[gene_key, 'accession']]
                if direction == 'normal': rel_start = rel_start * -1
                gene_key_hits[gene_key].append((queryfamkey, rel_start))

# now gene_key_hits is a dict that looks like this:
# {  404: [(5, -1703),  (2, -320),  (1, 1494)],
#    405: [(1, 1118),   (2, -1424), (5, -1855)],
#   ...
# }
# where each number corresponds to the fam's index in 'queryfams'
# note! The source enzyme (StzF) is NOT YET INCLUDED! it will be added next
# this next block goes thru each entry in GENE_KEY_HITS and
# adds the SSN reference hit (i.e. StzF-homolog in this case) (assigned '0')
# then it sorts the entries in GENE_KEY_HITS by 'rel_start', which should
# make all the BGCs displayed in the same order

for gene_key, hit in gene_key_hits.items():
    hit.append((0, organism_df.loc[gene_key, 'rel_start']))
    hit.sort(key=lambda x: x[1], reverse=True)
    hit = [x[0] for x in hit]
    if hit not in genotypes:
        genotypes.append(hit)
        genotype_count[genotypes.index(hit)] = 0
        genotype_genes[genotypes.index(hit)] = []
    genotype_count[genotypes.index(hit)] += 1
    genotype_genes[genotypes.index(hit)].append(gene_key)
    pass

# now you have genotype_count which is a dictionary containing all the
# different genotypes found, and the number of times they appear.
# genotype_genes is another dictionary which tells you the gene keys that match
# each of those genotypes.

# next: get rid of genotypes that only show up once:

genotype_count = {k: v for k, v in genotype_count.items() if v > 2}
genotype_genes = {k: v for k, v in genotype_genes.items() if k in genotype_count}

# next: generate a .csv containing the id and genotype # for each organism!
# (this is for importing into cytoscape)
for genotype, gene_keys in genotype_genes.items():
    for gene_key in gene_keys:
        cytotable_genotypes[gene_key].append(genotype)
    pass

# this next block prepares the table for export to cytoscape..
cytotable = pd.DataFrame.from_dict(cytotable_genotypes,
                                   orient='index',
                                   columns=['accession', 'genotype'])

cytotable.to_csv(path_or_buf='cytoscape_output.csv')

# this generates an html file with a graphical report. See legendhtmlcompile.py
legendhtmlcompile.export_data_report('legend.html', queryfams, queryfam_labels,
                                     genotypes, genotype_count)

# TODO: - import the 'queryfams' from a .csv provided by the user
#       - also format should be like this:
#       [label] [pfam synonym1] [pfam synonym2] [pfam synonym3] ..
#       - export cytotable (df), and export a table that translates the
#           genotypes into pfam labels. Maybe later generate a graphic instead?
