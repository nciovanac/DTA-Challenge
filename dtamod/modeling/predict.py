import torch
from torch import nn
from torch import LongTensor, cat, tensor, is_tensor
import torch.nn.functional as F
import torch.autograd as autograd
from torch.utils.data import Dataset, DataLoader
import argparse


import polars as pl

from train_early_stop import DeepDTA, DTADatasetLoader, lt, st, proteins, make_context_vector, pad_str


device = torch.device("cpu") # For inference this was a little faster than GPU
print(device)


def runInference(filepath):
    """
    Loads a trained model from state dictionary and predicts the output on a supplied input file.
    The load and prepare scripts must be triggered before this script
    """
    

    # Neural network and data loaders defined in the train_early_stop.py file 
    deepdta = DeepDTA()
    deepdta.load_state_dict(torch.load("../../models/model_weights.pt",map_location = device, weights_only=True)) # In case you need to do inference on CPU 
    deepdta.eval()

    test_dataset = DTADatasetLoader("../../data/processed/TEST_input.csv")
    test_dataloader = DataLoader(test_dataset, batch_size = 256, shuffle = False)


    # Loop over batches of test dataloader to generate predicted KIBA score for all drug-targt pairs
    output = []
    with torch.no_grad():
        for data in test_dataloader:
            smiles, prots, kba = data.values()
            pbinding = deepdta(smiles.to(device), prots.to(device))
            output.extend(pbinding.squeeze(-1).tolist())



    # Adds predictions and metrics to test dataframe
    df_out = pl.read_csv(filepath)

    df_out = df_out.with_columns(pl.Series(name="predicted_kiba_score", values = output))

    df_out = df_out.with_columns(
        (abs((pl.col("predicted_kiba_score") - pl.col("kiba_score")))).alias("abs_error")

    )

    df_out = df_out.with_columns(
        (pl.col("abs_error")**2).alias("sq_error")
        )


    df_out.write_csv("../../data/processed/TEST_output.csv")

    print(df_out.select(pl.mean("abs_error","sq_error")))

#    print(pbinding.squeeze(-1).tolist())




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Loads a trained model from state dictionary and predicts the output on a supplied input file. Predictions and error metrics will be added as additional columns and saved to a separate file. MAE/MSE on the test set will be printed to stout. The load and prepare scripts must be triggered before this script.")
    parser.add_argument("filename", help = "Path to test file")

    args = parser.parse_args()
    
    runInference(args.filename)
