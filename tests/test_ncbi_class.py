import unittest
import pandas as pd

from ncbi import NCBIdata

class TestNCBI(unittest.TestCase):
    def setUp(self) -> None:
        self.ncbi = NCBIdata()
    
    def tearDown(self) -> None:
        return super().tearDown()

    def test_data_loading(self):
        
        self.ncbi.load()

        # df is not None
        self.assertIsNotNone(self.ncbi.genome_df)
        self.assertTrue(self.ncbi.genome_df.shape[0] > 0) # number of rows
        self.assertTrue(self.ncbi.genome_df.shape[1] > 0) # number of columns
        self.assertTrue(self.ncbi.genome_df.shape[1] == len(self.ncbi.column_names)-1)

    
    


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