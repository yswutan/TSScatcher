import torch
from torch.utils.data import DataLoader, Dataset
import pandas as pd
from torch.nn.utils.rnn import pad_sequence

class CustomDataset(Dataset):
    def __init__(self, sequences):
        self.sequences = sequences

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        return self.sequences[idx]


def one_hot_encode(sequences):
    encoding = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1],
                "N": [0.25, 0.25, 0.25, 0.25],
                'a': [1, 0, 0, 0], 'c': [0, 1, 0, 0], 'g': [0, 0, 1, 0], 't': [0, 0, 0, 1],
                "n": [0.25, 0.25, 0.25, 0.25]}
    encoded_sequences = []
    for sequence in sequences:
        encoded_seq = [encoding[base] for base in sequence]
        encoded_sequences.append(encoded_seq)
    return encoded_sequences


def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = []
        for line in file:
            if line.startswith('>'):
                continue
            sequences.append(line.strip())
    df = pd.DataFrame(sequences, columns=["Sequence"])
    print(df)
    return df


def custom_collate_fn(batch):
    batch = [torch.tensor(item, dtype=torch.float32) for item in batch]

    padded_batch = pad_sequence(batch, batch_first=True)

    return padded_batch


def process_data_eval(data_path, batch_size):
    df = read_fasta(data_path)
    temp_seq = list(df["Sequence"].values)
    encoded_sequences = one_hot_encode(temp_seq)
    X = torch.tensor(encoded_sequences, dtype=torch.float32)

    dataset =CustomDataset(X)
    evl_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    return evl_loader
