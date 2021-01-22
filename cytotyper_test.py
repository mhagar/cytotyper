import pandas as pd
import sqlite3
import os

gap_unit = 300  # Number of bps before adding a 'gap' in the BGC legend
min_occurence = 1  # Number of times a BGC type has to occur or it'llget tossed

# See "datastructures.pptx" for a diagram of the datastructures.


def importdata(rawdata):
    # connect to sqlite file
    conn = sqlite3.connect(os.path.join('uploads', rawdata.filename))

    sqlquery = "id, num, family, ipro_family, rel_start, rel_stop, \
                direction, gene_key, start"
    main_df = pd.read_sql_query(f"SELECT {sqlquery} FROM neighbors", conn)

    # 'attributes' table in db, index corresponds to gene key in 'neighbors'!
    # i.e. gene_key 573, which refers to kutzneria bgc, the stzF homolog is
    # found in 'attributes' row 573.

    sqlquery = "sort_key, accession, num, family,\
                 ipro_family, direction, start, strain"
    organism_df = pd.read_sql_query(f"SELECT {sqlquery} FROM attributes", conn,
                                    index_col='sort_key')

    # 'AHMAD-BORSCHT-CAR' -> ['AHMAD', 'BORSCHT', 'CAR']
    def split_dash(x):
        x = x.split('-')
        return x

    main_df['family'] = main_df['family'].apply(split_dash)
    main_df['ipro_family'] = main_df['ipro_family'].apply(split_dash)

    # combine 'family' and 'ipro_family' columns
    main_df['family'] += main_df['ipro_family']
    main_df.drop(labels='ipro_family', axis=1, inplace=True)
    return main_df, organism_df


def gene_query():
    # Just a test query for now.
    global query
    query = {1: ['PF06325', 'IPR025799', 'IPR029063'],
             2: ['PF00355', 'IPR036922', 'IPR017941'],
             3: ['PF13535', 'IPR011761'],
             4: ['PF01425', 'IPR036928', 'IPR023631']}

    global query_legend
    query_legend = {1: 'PrmA',
                    2: 'Reiske',
                    3: 'ATP-grasp',
                    4: 'Amidase'}


class gene_hit:
    def __init__(self, pfamkey, direction,
                 rel_start, rel_stop):
        self.pfamkey = pfamkey
        self.direction = direction
        self.rel_start = rel_start
        self.rel_stop = rel_stop
        self.length = abs(rel_start - rel_stop)


class gap:
    def __init__(self, length):
        self.length = length


def groupdata(df, gene_key):
    # 'grouped_df.get_group(573)' returns a df with only the rows with
    #  genekey = 573
    if gene_key in df['gene_key'].values:
        grouped_df = df.groupby(by="gene_key")
        return grouped_df.get_group(gene_key).reset_index()
    else:
        return 0


class cluster:
    def __init__(self, data, gene_key, accession, organism_id):
        self.data = data
        self.gene_key = gene_key
        self.accession = accession
        self.organism_id = organism_id
        self.summary = []
        self.hash = []


def parse_cluster(cluster, attributes):
    # for every row in the gene cluster table..
    data = cluster.data
    if isinstance(data, pd.DataFrame) is False:
        return "Null"

    # Generating cluster.summary:
    for i in range(0, len(data)):
        hit_flag = 0

        # if you just passed the row where "stzF" would be..
        if data.loc[i, 'start'] > attributes['start']:
            # AND if you didnt already place "stzF" in this cluster.summary..
            if 0 not in [x.pfamkey for x in cluster.summary
                         if isinstance(x, gene_hit)]:
                # add an entry for "stzF".
                hit = gene_hit(pfamkey=0,
                               direction='normal',
                               rel_start=0,
                               rel_stop=0)
                cluster.summary.append(hit)

        # for every item in the query .
        for pfamkey, pfams in query.items():
            # if there's a match between a label's pfams and the cluster table-
            if [x for x in pfams if x in data.loc[i, 'family']] != []:
                # create a 'gene_hit' object..
                hit = gene_hit(pfamkey=pfamkey,
                               direction=data.loc[i, 'direction'],
                               rel_start=data.loc[i, 'rel_start'],
                               rel_stop=data.loc[i, 'rel_stop'])

                cluster.summary.append(hit)
                hit_flag = 1

        # if this row didn't match anything in the query..
        if hit_flag == 0:
            gap_length = abs(data.loc[i, 'rel_start']
                             - data.loc[i, 'rel_stop'])
            # if there's already a gap object at the end of clusterlist..
            if len(cluster.summary) > 0 \
                    and isinstance(cluster.summary[-1], gap):
                # extend its length property.
                cluster.summary[-1].length += gap_length
            else:
                # add a gap object.
                cluster.summary.append(gap(gap_length))

    # Generating cluster.hash:
    for k, item in enumerate(cluster.summary):
        if isinstance(item, gap):
            # If you're at the beginning OR end of the cluster.summary..
            if k == 0 or k+1 == len(cluster.summary):
                # Append ONE 'G'
                cluster.hash.append('G')
            else:
                # Otherwise, append a 'G' for every 'gap_unit' basepairs
                # default: 300
                for n in range(0, (item.length // gap_unit)):
                    cluster.hash.append('G')

        if isinstance(item, gene_hit):
            direction = 1
            if item.direction == 'complement':
                direction = -1
            cluster.hash.append(str(item.pfamkey * direction))


def parse_dataset(neighbor_table, attributes_table):
    parsed_dataset = {}
    data = attributes_table
    for i in range(1, len(data)+1):
        c = cluster(data=groupdata(neighbor_table, i),
                    gene_key=i,
                    accession=data['accession'].loc[i],
                    organism_id=data['strain'].loc[i])
        parse_cluster(c, data.loc[i])
        parsed_dataset[i] = c

    return parsed_dataset


class genotype:
    def __init__(self, string_hash):
        self.count = 0
        self.gene_keys = []
        self.accessions = []
        self.organisms = []
        self.string_hash = ""


def bin_hashes(dataset):

    # Go through all the clusters and create a new dict entry for every new
    # hash encountered. If it's already in the dict, increment the count by one
    unique_string_hashes = {}
    genotypes = {}

    for key, cluster in dataset.items():
        string_hash = ",".join(cluster.hash)

        # if this string-hash hasn't been encountered before..
        if string_hash not in unique_string_hashes:
            # create a new entry
            unique_string_hashes[string_hash] = 0
            genotypes[string_hash] = genotype(string_hash=string_hash)

        unique_string_hashes[string_hash] += 1
        genotypes[string_hash].count += 1
        genotypes[string_hash].gene_keys.append(cluster.gene_key)
        genotypes[string_hash].accessions.append(cluster.accession)
        genotypes[string_hash].organisms.append(cluster.organism_id)

    # Now go through the genotype dict and remove entries that only show up
    # once
    filtered_genotypes = {}
    for string_hash, item in genotypes.items():
        if item.count > min_occurence:
            filtered_genotypes[string_hash] = item

    # clear the <genotypes> dict
    genotypes = {}

    # replace it with a numbered version
    genotypes = dict(zip(range(0, len(filtered_genotypes)),
                         filtered_genotypes.values()))

    return genotypes


def generate_cytotable(genotypes):
    d = {}
    for key, item in genotypes.items():
        for accession in item.accessions:
            d[accession] = key

    return pd.DataFrame.from_dict(d, orient='index')


def analyze(filefield):
    """
    Input: a FileField object with a .filename property (string) that says
    what the file is called in the /uploads directory

    Returns: a tuple (genotypes dict, cytotable Dataframe)
    """
    gene_query()  # TODO: add a system for specifying query
    dataset = parse_dataset(importdata(filefield)[0],
                            importdata(filefield)[1])
    genotypes = bin_hashes(dataset)
    cytotable = generate_cytotable(genotypes)
    return genotypes, cytotable


# ### TEST GROUNDS
# ###


class test_input:
    def __init__(self, filename):
        self.filename = filename


tinput = test_input('StzF-truncated.sqlite')

analyze(tinput)

import pdb; pdb.set_trace()
#
# tinput = test_input('StzF-SSN-Nov-26-2020.sqlite')
# main_df = importdata(tinput)[0]
# organism_df = importdata(tinput)[1]
# gene_query()
# # testkutz = cluster(groupdata(main_df, 533).reset_index())
#
# import pdb; pdb.set_trace()
#
# test = parse_dataset(main_df, organism_df)
#
# import pdb; pdb.set_trace()
#
# genotypes = bin_hashes(test)
#
# import pdb; pdb.set_trace()

# ###TODO: export .csv for cytoscape
# ###TODO: generate .html report of said genotypes above.
