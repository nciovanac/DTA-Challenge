import torch
from torch import nn
from torch import LongTensor, cat, tensor, is_tensor
import torch.nn.functional as F
import torch.autograd as autograd
from torch.utils.data import Dataset, DataLoader
import torch.optim as optim

import polars as pl


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)



class DeepDTA(nn.Module):

    def __init__(self):
        super().__init__()

        # Conv1D, filters = numfilers, kernel_size = filter_length1, activation = relu, padding = "valid", strides = 1)

        nfilter = 32
        filter_length1 = 4
        filter_length2 = 6
        self.smiles_embedding = nn.Embedding(40, 100) # num chars, max_len 
        self.smiles_encode1 = nn.Conv1d(100, nfilter, filter_length1)
        self.smiles_encode2 = nn.Conv1d(nfilter, nfilter*2, filter_length1)
        self.smiles_encode3 = nn.Conv1d(nfilter*2, nfilter*3, filter_length1)
        self.smiles_latent = nn.AdaptiveMaxPool1d(1)

        self.protein_embedding = nn.Embedding(26, 1000)
        self.protein_encode1 = nn.Conv1d(1000, nfilter, filter_length2)
        self.protein_encode2 = nn.Conv1d(nfilter, nfilter*2, filter_length2)
        self.protein_encode3 = nn.Conv1d(nfilter*2, nfilter*3, filter_length2)
        self.protein_latent = nn.AdaptiveMaxPool1d(1)

        self.fc1 = nn.Linear(192,1024)
        self.fc2 = nn.Linear(1024,1024)
        self.fc3 = nn.Linear(1024,512)
        self.fc4 = nn.Linear(512, 1)

        self.dropout = nn.Dropout()

    def forward(self, x, y):

        x = F.relu(self.smiles_embedding(x))

        x = F.relu(self.smiles_encode1(x))

        x = F.relu(self.smiles_encode2(x))

        x = F.relu(self.smiles_encode3(x))

        x = F.relu(self.smiles_latent(x))


        y = F.relu(self.protein_embedding(y))

        y = F.relu(self.protein_encode1(y))
        
        y = F.relu(self.protein_encode2(y))
        
        y = F.relu(self.protein_encode3(y))
        
        y = F.relu(self.protein_latent(y))

        z = cat((x,y), dim = 1).squeeze(2)

        #z = cat((x,x), dim = 1).squeeze(2)

        z = F.relu(self.fc1(z))
        z = self.dropout(z)
        z = F.relu(self.fc2(z))
        z = self.dropout(z)
        z = F.relu(self.fc3(z))
        z = self.fc4(z)
        return z


#deepdta = DeepDTA()




#loss_fun = nn.MSELoss()

#optimizer = optim.Adam(deepdta.parameters(), lr = 0.001)




# Based on word2vec

# Dictionary of valid SMILES characters


lt = {
    "Br":"X",
    "Ni": "Y",
    "Cl" : "Z"}


st = ['K', 'b','(', '+', '9', 'S', 'F', '8', '.', '1', '5', '-', 'N', '/', '[', 'n', 'C', '7', ')', 'O', '%', '0', '#', '6', 'B', '\\', '@', ']', '4', '=', 'P', '2', '3', 's', 'I','H', "X","Y","Z", " "]

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


ex_prot = "MPEEVQTQDQPMETFAVQTFAFQAEIAQLMSLIYESLTDPSKLDSGK"
ex_smile = "CCN1C=NC2=C(N=CN=C21)N"
ex_binding = 10.0

def make_context_vector(context, word_to_ix):
    idxs = [word_to_ix[w] for w in context]
    tensor = torch.tensor(idxs, device = device)
    return autograd.Variable(tensor)

def pad_str(s, N):
    return(s.ljust(N)[:N])



smiles_vocab = set(st)
smiles_vocab_size = len(st)

prot_vocab = set(proteins)
prot_vocab_size = len(proteins)

smile_to_ix = {smile: i for i, smile in enumerate(smiles_vocab)}
prot_to_ix = {a: i for i, a in enumerate(prot_vocab)}
#data = []

#data.append((pad_str(ex_smile,50 ),pad_str(ex_prot,100), ex_binding))

# for epoch in range(5):
#     total_loss = 0
#     for smile, protein, binding in data:
#         smile_vec = make_context_vector(smile, smile_to_ix)
#         prot_vec = make_context_vector(protein, prot_to_ix)
#         smile_vec = smile_vec[None, :]
#         prot_vec = prot_vec[None,:]
        
        
#         deepdta.zero_grad()
#         pbinding = deepdta(smile_vec, prot_vec)
        
#         loss = loss_fun(tensor(binding).reshape([1,1]), pbinding)
#         loss.backward()
#         optimizer.step()
#         total_loss += loss.data
#         print(total_loss)






class DTADatasetLoader(Dataset):
    def __init__(self, filepath):
        self.dat = pl.read_csv(filepath)
        print(len(self.dat))

    def __len__(self):
        return len(self.dat)

    def __getitem__(self, idx):
        if is_tensor(idx):
            idx = id.tolist()
        
        smiles = pad_str(self.dat[idx, "smiles"], 100)
        prots = pad_str(self.dat[idx, "Sequence"], 1000)
        kba = self.dat[idx, "kiba_score"]


        smiles = make_context_vector(smiles, smile_to_ix)
        prots = make_context_vector(prots, prot_to_ix)
        kba = tensor(kba, device = device).float()

        sample = {
            "smiles" : smiles,
            "prots": prots,
            "kba" : kba
            }

        return sample



def main():

    dta_dataset = DTADatasetLoader("../../data/processed/TRAINING.csv")

    dataloader = DataLoader(dta_dataset, batch_size = 256, shuffle = True)


    train_dat = next(iter(dataloader))



    print("training model...")

    num_epochs = 100


    deepdta = DeepDTA()
    deepdta.to(device)

    loss_fun = nn.MSELoss()

    optimizer = optim.Adam(deepdta.parameters(), lr = 0.001)


    train_dat = DataLoader(dta_dataset, batch_size = 512, shuffle = True)

    for epoch in range(num_epochs):
        print(f'Epoch {epoch}/{num_epochs - 1}')
        print('-' * 10)
        for i, data in enumerate(train_dat):
            print("On batch {}".format(i))
            optimizer.zero_grad()
            smiles,prots,kba = data.values()

        

        
            pbinding = deepdta(smiles, prots)
        
            loss = loss_fun(kba, pbinding.squeeze(1))
            loss.backward()
            optimizer.step()

            #        print(float(loss.mean()))
            #        print(total_loss)




    print("finished!")

    # Save model to file

    torch.save(deepdta.state_dict(), "../../models/working_model2.pt")


if __name__ == "__main__":
    main()
