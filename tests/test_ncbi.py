import unittest
import pandas as pd

class TestNCBI(unittest.TestCase):
    def setUp(self) -> None:
        DATA_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
        self.column_names = ['#Organism/Name', 'TaxID', 'Size (Mb)', 'GC%', 'Replicons', 'Genes', 'Proteins', 'Release Date', 'Status', 'FTP Path']
        self.df = pd.read_table(DATA_URL, usecols=self.column_names)  # the column 15: Status
      
    
    def tearDown(self) -> None:
        return super().tearDown()

    def test_data_loading(self):
        
        # df is not None
        self.assertIsNotNone(self.df)
        self.assertTrue(self.df.shape[0] > 0) # number of rows
        self.assertTrue(self.df.shape[1] > 0) # number of columns
        self.assertTrue(self.df.shape[1] == len(self.column_names))

    
    def test_data_loading2(self):
        
        complete_genomes = self.df.loc[self.df['Status'] == "Complete Genome"]    # Select only complete genomes
        genome_df = complete_genomes.reset_index(drop = True) # reset index from 0
        
        # 'Status' == "Complete Genome" ???
        
        genome_df = genome_df.drop(['Status'], axis = 1) # remove the column "Status"
        genome_df = genome_df.rename(columns={"#Organism/Name": "Genome Name"}) # change the first column name
    
        # which testing?

    


    # def test_load_data(self):
        
    #     ## load genome data from NCBI 
    #     DATA_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
    #     df = pd.read_table(DATA_URL)  # the column 15: Status
    #     self.assertIsNotNone(df)

    #     column_names = ['#Organism/Name', 'TaxID', 'Size (Mb)', 'GC%', 'Replicons', 'Genes', 'Proteins', 'Release Date', 'Status', 'FTP Path']
    #     print(list(df.columns.values))
    #     print(df.shape)

    #     self.assertTrue(df.shape[0]>0)



    

    

if __name__ == '__main__':
    unittest.main()