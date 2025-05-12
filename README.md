# fastBJM new new

Things in the dat object
- age_intervals: definition of age intervals (as defined by user)
- nGs: number of age intervals
- country: either country names or country ids (not necessary consecutively though e.g. can be 11, 12, 14 etc)
- nctys: number of countries included in the study
- cty_id: country ID from 1 to nctys (i.e. as.numeric(as.factor(country)))

- jkg_index: a string transition-age index in the format of '12g1' (this is also how the baseline hazards are named in the params object)
- jkc_index: a string transition-country index in the format of 'cty1_12'
- jkgc_index: a string transition-age_country index in the format of 'cty1_12g1'

- X_spec (to add): an nXs-by-2 array with the names of the X variables (column 1) and their types ('cnt','cat') where cnt=continuous and cat=binary/categorical; the variable names MUST be the same as those in the data file
- X: a matrix of fixed effect values; users need to provide the column names (e.g. 'edu1','edu2','edu3' for three education levels)
- allowable_transitions: permissible transitions as in the multistate model
- nats (to add): number of allowable transitions (i.e. length(allowable_transitions))
- jk_index (to add): a string transition index in the form of '12'
- trans_indicators: a dummy version of jk_index with column names being allowable_transitions

- nQs: number of quadrature points (default=15 with min=1)
- Q
- 