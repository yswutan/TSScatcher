import argparse
import torch
from model import Net
from train_test import val_model
from utils_test import process_data_eval
import os
def save_indices_to_txt(indices, filename):
    try:
        with open(filename, 'w') as file:
            for index in indices:
                file.write(f"{index}\n")
    except IOError as e:
        print(f"Error writing to file {filename}: {e}")

def predict(data_path, batch, load_path):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    eval_loader = process_data_eval(data_path, batch)
    model = Net()
    model.to(device)

    try:
        model.load_state_dict(torch.load(load_path,map_location=torch.device("cpu"), weights_only=True),strict=True)
    except IOError as e:
        print(f"Error loading model from {load_path}: {e}")
        return None

    return val_model(model, eval_loader, device)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate model predictions on new data.")
    parser.add_argument("--input_file", required=True, help="Path to the input file.")
    parser.add_argument("--output_file", required=True, help="Path to the output file where results will be saved.")

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(current_dir, "tsscatcher.pth")
    batch_size = 256

    #print(f"Evaluating data from: {input_file}")
    #print(f"Saving results to: {output_file}")

    result = predict(input_file, batch_size, model_path)

    if result is not None:
        predict_probability = result[0]
        indices = [value.item() for value in predict_probability]
        save_indices_to_txt(indices, output_file)
        print(f"Deep learning model: Done.")
    else:
        print("Prediction failed.")
