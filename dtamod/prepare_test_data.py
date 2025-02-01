import polars as pl
import argparse


def convertLongTokens(smiles):
    """Attempts to convert the most common two character smiles into single character tokens.
    If any more exotic characters are found, e.g. W, they will be converted into a blank token
    """
    
    long_tokens = {
        "Br":"X",
        "Ni": "Y",
        "Cl" : "Z"}


    short_tokens = ['K', 'b','(', '+', '9', 'S', 'F', '8', '.', '1', '5', '-', 'N', '/', '[', 'n', 'C', '7', ')', 'O', '%', '0', '#', '6', 'B', '\\', '@', ']', '4', '=', 'P', '2', '3', 's', 'I','H', "X","Y","Z", " "]

    set_short_tokens = set(short_tokens)
    # Dictionary of prots


    proteins = [ "A",
             "C",
             "B",
             "E",
             "D",
             "G",
             "F",
             "I",
             "H",
             "K",
             "M",
             "L",
	     "O",
             "N",
             "Q",
             "P",
             "S",
             "R",
	     "U",
             "T",
             "W",             
	     "V",
             "Y",
             "X",            
	     "Z",
             " "]

    
    for lt, st in long_tokens.items():
        smiles = smiles.replace(lt, st)
    irregulars = set(smiles) - set_short_tokens
    if irregulars:
        for token in list(irregulars):
            smiles = smiles.replace(token, " ")
    return smiles
        



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Joins SMILES and protein sequence columns with dataset based on UniProt ID and Pubchem CID. You must run the 'load_test_data.py' script first in order to generate tables with ID -> SMILES/Sequence mappings. Processed files will be written to predetermined locations for ease of testing")
    parser.add_argument("filename", help = "Absolute path to test file")

    args = parser.parse_args()

    dat = pl.read_csv(args.filename)

    dat = dat.with_columns(pl.col("pubchem_cid").cast(pl.Int64, strict = False))

    # Load smiles mapping

    smiles_mapping = pl.read_csv("../data/interim/TEST_smiles_mapping.txt", separator = "\t", has_header = False, new_columns = ["pubchem_cid", "smiles"])

    # Join on dataset

    dat = dat.join(smiles_mapping, on = "pubchem_cid")

    dat = dat.with_columns(pl.col("smiles").map_elements(convertLongTokens, return_dtype = pl.String))

    # Load protein mapping

    protein_mapping = pl.read_csv("../data/processed/all_prots.csv")

    # JOin on dataset

    dat = dat.join(protein_mapping, on = "UniProt_ID")

    # Write to file

    test_location = "../data/processed/TEST_input.csv"
    dat.write_csv(test_location)

    print("Test data prepared and written to {}".format(test_location))

    quit()







