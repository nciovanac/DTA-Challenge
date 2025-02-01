# SFL Drug-Target Affinity Challenge

**Prepared for SFL Scientific Sr. Data Scientist Challenge**

Testing Instructions

Venv requirements have been exported to "requirements.txt". No GPU is required to run the inference pipeline. Place the test file in the "./data/raw" directory. From within the "dtamode" directory, load the SMILES and protein sequences with "python load_test_data.py <path_to_your_test_file>". Once the data is pulled, prepare the data by running "python prepare_test_data.py <path_to_your_test_file>". Finally, enter the "modeling" directory and run "python predict.py <path_to_your_test_file>. MAE and MSE on the test set will be printed to screen once available
