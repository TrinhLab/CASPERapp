import sys
from CSPRparser import CSPRparser
import pandas as pd


"""
This function takes one input currently: a CSPR file passed after the executable.
To run this program, simply execute the command:

    python3 main.py <path to CSPR file>

"""
class GoldenGateGuide():
    def __init__(self,cspr_path):
        try:
            self.cspr_path = cspr_path # Set path to CSPR file
            self.parser_obj = CSPRparser(cspr_path) # Create CSPR parser object
            self.parser_obj.read_first_lines() # Read basic info from CSPR file
            self.pos_tuples = self.get_pos_tuples(self.parser_obj.karystatsList) # Gets the query arguments for CSPRparser
            self.guide_list = [] # Initialize list that will hold all the gRNAs in the CSPR file
            for pos_tuple in self.pos_tuples: # Iterate through each chromosome and pull all the guides for it
                self.guide_list.append(self.parser_obj.read_targets(pos_tuple=pos_tuple, endo="spCas9")) # Append guides to guide_list
            self.IIS_enzymes = self.get_enzymes()


        except Exception as e:
            print(e)

    def get_pos_tuples(self, karyStats):
        pos_tuples = []
        for i, length in enumerate(karyStats):
            pos_tuples.append((i+1,1,length))
        return pos_tuples
    
    def get_enzymes(self):
        df = pd.read_excel("NEB_IIS_Enzymes.xlsx",header=0,engine="openpyxl") # Read IIS table
        df = df[df["Overhang Length"]>=3] # Filter out enzymes with overhangs less than 2 nt
        df.replace(u'\xa0',u'', regex=True, inplace=True) # Remove unicode space characters
        enzyme_dict = df.set_index('Enzyme')[['Recognition Sequence','Recognition Sequence Length','Overhang Length']].to_dict(orient='index')
        return enzyme_dict




""" Execute the code """
if __name__ == "__main__":
    main_obj = GoldenGateGuide(sys.argv[1])
