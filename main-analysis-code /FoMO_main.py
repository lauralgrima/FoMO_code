import load_all_data as lad

# for information about which experiments mice did etc. see the 'data_summary' word doc uploaded in the github repo

df_folder = '/Users/grimal/Dropbox/100hours/FoMO_DA/data/dataframes'

subject_IDs = ['6PG5','6PG6','6PG8','6PG9','6PG10','6PG11','6PG12','6PG15','6PG24','6PG25','6PG27','6PG28','6PG30','6PG31','6PG33','6PG35','6PO3','6PO4','6PO5',
               '6PO7','6PO8','6PO9','6PO10','dProb1','dProb3','gProb1','gProb2','ML11','ML12','ML13','ML14','ML15','ML17']
regions     = ['NAc','DMS']
tasks       = ['conc','prob']

data_dict = lad.load_data(df_folder, subject_IDs, tasks)
