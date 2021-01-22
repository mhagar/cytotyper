# parsed_dataset (dict)
  ## 1: <cluster>
  ### data (pandas.DataFrame)
  ### gene_key
  ### accession
  ### organism_id
  ### summary (list)
  #### gap
   ##### length (int)
  #### gene_hit
   ##### pfamkey
   ##### direction
   ##### rel_start
   ##### rel_stop
   ##### length
  #### gap
  #### ..
  #### gap
  ### hash (list)
  ## 2: <cluster>
  ## 3: <cluster>
  ## ..
  ## n: <cluster>

# genotypes (dict)
## 1: <genotype>
  ### count int
  ### gene_keys [list]
  ### accessions [list]
  ### organisms [list]
  ### string_hash "str"
## 2: <genotype>
## 3: <genotype>
## ..: <genotype>
## p: <genotype>
