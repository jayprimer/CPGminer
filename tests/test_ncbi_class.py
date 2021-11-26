import unittest
import pandas as pd
from app_v2_1_1_final import load_data
import pprint

from ncbi import NCBIdata

class TestNCBI(unittest.TestCase):
    def setUp(self) -> None:
        self.ncbi = NCBIdata()
    
    def tearDown(self) -> None:
        return super().tearDown()

    def load_and_print(self, filename):
        df = self.ncbi.load_df(filename)
        print(f'{filename}: {df.shape}, {df.columns}')
        print(df.tail())

    def test_data_loading(self):
        
        self.ncbi.load()
        # self.ncbi.load_from_ncbi()

        # self.ncbi.processing()

        # self.load_and_print('original')
        # self.load_and_print('step1')
        # self.load_and_print('step2')
        # self.load_and_print('step3')
        # self.load_and_print('step4')
        # self.load_and_print('step5')
        
        pprint.pprint(self.ncbi.tax_items.keys())


        # print(self.ncbi.genome_df.shape)
        # print(self.ncbi.genome_df.columns)
        # print(self.ncbi.genome_df.tail())

        # df is not None
        self.assertIsNotNone(self.ncbi.genome_df)
        self.assertTrue(self.ncbi.genome_df.shape[0] > 0) # number of rows
        self.assertTrue(self.ncbi.genome_df.shape[1] > 0) # number of columns
        # self.assertTrue(self.ncbi.genome_df.shape[1] == len(self.ncbi.column_names)-1)

    
    


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