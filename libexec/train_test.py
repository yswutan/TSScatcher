from torch.utils.data import DataLoader
import torch.nn.functional as F
from tqdm.auto import tqdm
import torch
"""
train and test function part
"""

class EarlyStopping:
    def __init__(self, patience=5, delta=0):
        self.patience = patience
        self.delta = delta
        self.counter = 0
        self.best_loss = float('inf')
        self.early_stop = False

    def __call__(self, loss):
        if loss < self.best_loss - self.delta:
            self.best_loss = loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True


def save_model(model, path):
    torch.save(model.state_dict(), path)


def accuracy_fn(y_true, y_pred):
    correct = torch.eq(y_true, y_pred).sum().item()
    acc = (correct / len(y_pred)) * 100
    return acc


def train(model, train_dataloader, epoch, optimizer, loss_fn, save_model_path, device='cuda'):
    epochs = epoch
    early_stopping = EarlyStopping(patience=5, delta=0.01)
    for epoch in tqdm(range(epochs)):
        print(f"Epoch:{epoch}")
        train_loss = 0
        for batch, (X, y) in enumerate(train_dataloader):
            model.train()
            X, y = X.to(device), y.to(device)
            y_pred = model(X)
            loss = loss_fn(y_pred, y)
            train_loss += loss
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            if batch % 1000 == 0:
                print(f"Looked at {batch * len(X)}/{len(train_dataloader.dataset)} samples")
        train_loss /= len(train_dataloader)
        print(f'\nTrain loss:{train_loss:.5f}')
        early_stopping(train_loss)
        if early_stopping.early_stop:
            print("Early stopping!")
            save_model(model, save_model_path)
            break
    if not early_stopping.early_stop:
        save_model(model, save_model_path)


def eval_model(model: torch.nn.Module, data_loader: torch.utils.data.DataLoader, loss_fn: torch.nn.Module, accuracy_fn,
               device):
    """
    model
    data_loader
    loss_fn
    accuracy_fn
    return
    """
    true_labels = []
    predicted_probs = []
    predicted_label_list = []
    loss, acc = 0, 0
    model.eval()
    with torch.inference_mode():
        for X, y in data_loader:
            X, y = X.to(device), y.to(device)
            y_pred = model(X)
            loss_temp = loss_fn(y_pred, y)
            loss += loss_temp
            acc += accuracy_fn(y_true=y, y_pred=y_pred.argmax(dim=1))
            predicted_labels = y_pred.argmax(dim=1)
            predicted_label_list.extend(predicted_labels)
            predicted_probs.extend(F.softmax(y_pred, dim=1)[:, 1].cpu().numpy())  # Probability of positive class
            true_labels.extend(y.cpu().numpy())
        loss /= len(data_loader)
        acc /= len(data_loader)
    return {
        "model_name": model.__class__.__name__,
        "model_loss": loss.item(),
        "model_acc": acc,
        "model_predicted_proba": predicted_probs,
        "model_predicted_label_list":predicted_label_list,
        "model_true_labels": true_labels}


def val_model(model, val_data, device):
    model.eval()
    predicted_probs = []
    predicted_label_list = []
    with torch.inference_mode():
        for _,X in enumerate(val_data):
            X = X.to(device)
            y_pred = model(X)
            predicted_probs.extend(F.softmax(y_pred, dim=1)[:, 0].cpu().numpy())  # Probability of positive class
            predicted_labels = y_pred.argmax(dim=1)
            predicted_label_list.extend(predicted_labels)
    return predicted_probs, predicted_label_list
