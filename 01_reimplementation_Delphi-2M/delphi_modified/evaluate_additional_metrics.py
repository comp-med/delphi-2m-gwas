'''
This script is adapted from ./evaluate_auc.py
Changes:
- Replaces AUC-only evaluation with:
    precision, recall, f1-score, AUPRC, balanced accuracy
- Keeps disease chunking, sex stratification, and age-bracket evaluation
- Uses sigmoid(logits) >= 0.5 as the default threshold for hard predictions
- Optionally supports bootstrap averaging for the metrics
'''
import warnings
import torch
from model import DelphiConfig, Delphi
from tqdm import tqdm
import pandas as pd
import numpy as np
import argparse
from utils import get_batch, get_p2i
from pathlib import Path
from scipy.special import expit

from sklearn.metrics import (
    precision_score,
    recall_score,
    f1_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss
)

def get_common_diseases(delphi_labels, filter_min_total=100):
    chapters_of_interest = [
        "I. Infectious Diseases",
        "II. Neoplasms",
        "III. Blood & Immune Disorders",
        "IV. Metabolic Diseases",
        "V. Mental Disorders",
        "VI. Nervous System Diseases",
        "VII. Eye Diseases",
        "VIII. Ear Diseases",
        "IX. Circulatory Diseases",
        "X. Respiratory Diseases",
        "XI. Digestive Diseases",
        "XII. Skin Diseases",
        "XIII. Musculoskeletal Diseases",
        "XIV. Genitourinary Diseases",
        "XV. Pregnancy & Childbirth",
        "XVI. Perinatal Conditions",
        "XVII. Congenital Abnormalities",
        "Death",
    ]
    labels_df = delphi_labels[
        delphi_labels["ICD-10 Chapter (short)"].isin(chapters_of_interest) * (delphi_labels["count"] > filter_min_total)
    ]
    return labels_df["index"].tolist()

# convert logits to probabilities
def sigmoid_np(x):
    #print(f"logits: {x}")
    #print(f"probs: {expit(x)}")
    return expit(x)

# compute requested metrics instead of AUC
def compute_metrics(y_true, y_score, threshold=0.5):
    """
    Compute:
      - auprc
      - balanced accuracy
      - brier score

    Args:
        y_true: binary labels, 0/1
        y_score: predicted probabilities
        threshold: threshold for hard predictions
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)

    # Need both classes present for AUPRC / balanced metrics to be meaningful
    if len(np.unique(y_true)) < 2:
        return None

    y_pred = (y_score >= threshold).astype(int)

    return {
        #"precision": precision_score(y_true, y_pred),
        #"recall": recall_score(y_true, y_pred),
        #"f1": f1_score(y_true, y_pred),
        "auprc": average_precision_score(y_true, y_score),
        "balanced_acc": balanced_accuracy_score(y_true, y_pred),
        "brier": brier_score_loss(y_true, y_score)
    }

# bootstrap now resamples generic metrics, not AUC
def bootstrap_metrics(y_true, y_score, threshold=0.5, n_bootstrap=1000, seed=1337):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    out = []

    for bootstrap_idx in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        y_true_b = y_true[idx]
        y_score_b = y_score[idx]

        if len(np.unique(y_true_b)) < 2:
            continue

        metrics = compute_metrics(y_true_b, y_score_b, threshold=threshold)
        if metrics is None:
            continue

        metrics["bootstrap_idx"] = bootstrap_idx
        out.append(metrics)

    return out

def get_calibration_metrics(j, k, d, p, offset=365.25, age_groups=range(45, 80, 5), precomputed_idx=None, n_bootstrap=1, threshold=0.5, seed=1337):
    age_step = age_groups[1] - age_groups[0]

    # Indexes of cases with disease k
    wk = np.where(d[2] == k)

    if len(wk[0]) < 2:
        return None

    # For controls, we need to exclude cases with disease k
    wc = np.where((d[2] != k) * (~(d[2] == k).any(-1))[..., None])

    wall = (np.concatenate([wk[0], wc[0]]), np.concatenate([wk[1], wc[1]]))  # All cases and controls

    # We need to take into account the offset t and use the tokens for prediction that are at least t before the event
    if precomputed_idx is None:
        pred_idx = (d[1][wall[0]] <= d[3][wall].reshape(-1, 1) - offset).sum(1) - 1
    else:
        pred_idx = precomputed_idx[wall]  # It's actually much faster to precompute this

    z = d[1][(wall[0], pred_idx)]  # Times of the tokens for prediction
    z = z[pred_idx != -1]

    zk = d[3][wall]  # Target times
    zk = zk[pred_idx != -1]

    # use logits -> probabilities instead of using raw score for AUC only
    logits = p[..., j][(wall[0], pred_idx)]
    logits = logits[pred_idx != -1]
    scores = sigmoid_np(logits)

    wk = (wk[0][pred_idx[: len(wk[0])] != -1], wk[1][pred_idx[: len(wk[0])] != -1])
    p_idx = wall[0][pred_idx != -1]

    out = []

    for i, aa in enumerate(age_groups):
        a = np.logical_and(z / 365.25 >= aa, z / 365.25 < aa + age_step)
        # Optionally, add extra filtering on the time difference, for example:
        # a *= (zk - z < 365.25)
        selected_groups = p_idx[a]
        perm = np.random.permutation(len(selected_groups))
        _, indices = np.unique(selected_groups[perm], return_index=True)
        indices = perm[indices]
        selected = np.zeros(np.sum(a), dtype=bool)
        selected[indices] = True
        a[a] = selected

        control = scores[len(wk[0]) :][a[len(wk[0]) :]]
        case = scores[: len(wk[0])][a[: len(wk[0])]]

        if len(control) == 0 or len(case) == 0:
            continue

        # CHANGED: build binary labels for classification metrics
        y_true = np.concatenate([
            np.ones(len(case), dtype=int),
            np.zeros(len(control), dtype=int),
        ])
        y_score = np.concatenate([case, control])

        if n_bootstrap > 1:
            metrics_bootstrapped = bootstrap_metrics(
                y_true,
                y_score,
                threshold=threshold,
                n_bootstrap=n_bootstrap,
                seed=seed,
            )
            for item in metrics_bootstrapped:
                out_item = {
                    "token": k,
                    "age": aa,
                    "n_healthy": len(control),
                    "n_diseased": len(case),
                }
                out_item.update(item)
                out.append(out_item)
        else:
            metrics = compute_metrics(y_true, y_score, threshold=threshold)
            if metrics is None:
                continue

            out_item = {
                "token": k,
                "age": aa,
                "n_healthy": len(control),
                "n_diseased": len(case),
            }
            out_item.update(metrics)
            out.append(out_item)
    return out


# New internal function that performs the metrics evaluation pipeline.
def evaluate_metrics_pipeline(
    model,
    d100k,
    output_path,
    delphi_labels,
    diseases_of_interest=None,
    filter_min_total=100,
    disease_chunk_size=200,
    age_groups=np.arange(40, 80, 5),
    offset=0.1,
    batch_size=128,
    device="cpu",
    seed=1337,
    n_bootstrap=1,
    threshold=0.5,
    meta_info={},
):
    """
    Runs the evaluation pipeline using:
      auprc, balanced_acc, brier

    Args:
        model (torch.nn.Module): The loaded model set to eval().
        d100k (tuple): Data batch from get_batch.
        delphi_labels (pd.DataFrame): DataFrame with label info (token names, etc. "delphi_labels_chapters_colours_icd.csv").
        output_path (str | None): Directory where CSV files will be written. If None, files will not be saved.
        diseases_of_interest (np.ndarray or list, optional): If provided, these disease indices are used.
        filter_min_total (int): Minimum total token count to include a token.
        disease_chunk_size (int): Maximum chunk size for processing diseases.
        age_groups (np.ndarray): Age groups to use in calibration.
        offset (float): Offset used in get_calibration_metrics.
        batch_size (int): Batch size for model forwarding.
        device (str): Device identifier.
        seed (int): Random seed for reproducibility.
        n_bootstrap (int): Number of bootstrap samples. (1 for no bootstrap)
        threshold (float): 0.5
    Returns:
        tuple: (df_metrics_unpooled, df_metrics, df_both) DataFrames.
    """

    assert n_bootstrap > 0, "n_bootstrap must be greater than 0"

    # Set random seeds
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)

    # Get common diseases
    if diseases_of_interest is None:
        diseases_of_interest = get_common_diseases(delphi_labels, filter_min_total)

    # Split diseases into chunks for processing
    num_chunks = (len(diseases_of_interest) + disease_chunk_size - 1) // disease_chunk_size
    diseases_chunks = np.array_split(diseases_of_interest, num_chunks)

    # Precompute prediction indices for calibration
    pred_idx_precompute = (d100k[1][:, :, np.newaxis] < d100k[3][:, np.newaxis, :] - offset).sum(1) - 1

    all_metrics = []
    tqdm_options = {"desc": "Processing disease chunks", "total": len(diseases_chunks)}
    for disease_chunk_idx, diseases_chunk in tqdm(enumerate(diseases_chunks), **tqdm_options):
        p100k = []
        model.to(device)
        with torch.no_grad():
            # Process the evaluation data in batches
            for dd in tqdm(
                zip(*[torch.split(x, batch_size) for x in d100k]),
                desc=f"Model inference, chunk {disease_chunk_idx}",
                total=d100k[0].shape[0] // batch_size + 1,
            ):
                dd = [x.to(device) for x in dd]
                outputs = model(*dd)[0].cpu().detach().numpy()
                # Keep only the columns corresponding to the current disease chunk
                p100k.append(outputs[:, :, diseases_chunk].astype("float16"))  # enough to store logits, but not rates
        p100k = np.vstack(p100k)

        # Loop over each disease (token) in the current chunk, sexes separately
        for sex, sex_idx in [("female", 2), ("male", 3)]:
            sex_mask = ((d100k[0] == sex_idx).sum(1) > 0).cpu().detach().numpy()
            p_sex = p100k[sex_mask]
            d100k_sex = [d_[sex_mask].cpu().detach().numpy() for d_ in d100k]
            precomputed_idx_subset = pred_idx_precompute[sex_mask].cpu().detach().numpy()
            for j, k in tqdm(
                list(enumerate(diseases_chunk)), desc=f"Processing diseases in chunk {disease_chunk_idx}, {sex}"
            ):
                # Get calibration metrics for the current disease token.
                out = get_calibration_metrics(
                    j,
                    k,
                    d100k_sex,
                    p_sex,
                    age_groups=age_groups,
                    offset=offset,
                    precomputed_idx=precomputed_idx_subset,
                    n_bootstrap=n_bootstrap,
                    threshold=threshold,
                    seed=seed,
                )
                if out is None:
                    # print(f"No data for disease {k} and sex {sex}")
                    continue
                for out_item in out:
                    out_item["sex"] = sex
                    all_metrics.append(out_item)

    df_metrics_unpooled = pd.DataFrame(all_metrics)

    for key, value in meta_info.items():
        df_metrics_unpooled[key] = value

    delphi_labels_subset = delphi_labels[['index', 'ICD-10 Chapter (short)', 'name', 'color', 'count']]
    df_metrics_unpooled_merged = df_metrics_unpooled.merge(delphi_labels_subset, left_on="token", right_on="index", how="inner")

    def aggregate_age_brackets(group):
        return pd.Series({
            #'precision': group['precision'].mean(),
            #'recall': group['recall'].mean(),
            #'f1': group['f1'].mean(),
            'auprc': group['auprc'].mean(),
            'balanced_acc': group['balanced_acc'].mean(),
            'brier': group['brier'].mean(),
            'n_samples': len(group),
            'n_diseased': group['n_diseased'].sum(),
            'n_healthy': group['n_healthy'].sum(),
        })

    print('Aggregating precision / recall / f1 / auprc / balanced_acc across age brackets..')
    
    df_metrics = df_metrics_unpooled.groupby(["token"]).apply(aggregate_age_brackets).reset_index()
    df_metrics_merged = df_metrics.merge(delphi_labels, left_on="token", right_on="index", how="inner")

    if output_path is not None:
        Path(output_path).mkdir(exist_ok=True, parents=True)
        # CHANGED: new output filenames
        df_metrics_merged.to_parquet(f"{output_path}/df_metrics.parquet", index=False)
        df_metrics_unpooled_merged.to_parquet(f"{output_path}/df_metrics_unpooled.parquet", index=False)

    return df_metrics_unpooled_merged, df_metrics_merged


def main():
    parser = argparse.ArgumentParser(description="Evaluate metrics")
    parser.add_argument("--input_path", type=str, help="Path to the dataset")
    parser.add_argument("--output_path", type=str, help="Path to the output")
    parser.add_argument("--model_ckpt_path", type=str, help="Path to the model weights")
    parser.add_argument("--no_event_token_rate", type=int, help="No event token rate")
    parser.add_argument(
        "--health_token_replacement_prob", default=0.0, type=float, help="Health token replacement probability"
    )
    parser.add_argument("--dataset_subset_size", type=int, default=-1, help="Dataset subset size for evaluation")
    parser.add_argument("--n_bootstrap", type=int, default=1, help="Number of bootstrap samples")
    # Optional filtering/chunking parameters:
    parser.add_argument("--filter_min_total", type=int, default=100, help="Minimum total count to filter tokens")
    parser.add_argument("--disease_chunk_size", type=int, default=200, help="Chunk size for processing diseases")
    # add threshold argument for precision/recall/f1/balanced accuracy
    parser.add_argument("--threshold", type=float, default=0.5, help="Threshold for hard predictions")

    args = parser.parse_args()

    input_path = args.input_path
    output_path = args.output_path
    no_event_token_rate = args.no_event_token_rate
    health_token_replacement_prob = args.health_token_replacement_prob
    dataset_subset_size = args.dataset_subset_size

    # Create output folder if it doesn't exist.
    Path(output_path).mkdir(exist_ok=True, parents=True)

    device = "cuda"
    seed = 1337

    # Load model checkpoint and initialize model.
    ckpt_path = args.model_ckpt_path
    checkpoint = torch.load(ckpt_path, map_location=device)
    conf = DelphiConfig(**checkpoint["model_args"])
    model = Delphi(conf)
    state_dict = checkpoint["model"]
    model.load_state_dict(state_dict)
    model.eval()
    model = model.to(device)

    # Load training and validation data.
    val = np.fromfile(f"{input_path}/val.bin", dtype=np.uint32).reshape(-1, 3).astype(np.int64)

    val_p2i = get_p2i(val)

    if dataset_subset_size == -1:
        dataset_subset_size = len(val_p2i)

    # Get a subset batch for evaluation.
    d100k = get_batch(
        range(dataset_subset_size),
        val,
        val_p2i,
        select="left",
        block_size=80,
        device=device,
        padding="random",
        no_event_token_rate=no_event_token_rate,
        health_token_replacement_prob=health_token_replacement_prob,
    )

    # Load labels (external) to be passed in.
    delphi_labels = pd.read_csv("delphi_labels_chapters_colours_icd.csv")

    # call the new metrics pipeline instead of evaluate_metrics_pipeline
    df_metrics_unpooled, df_metrics_merged = evaluate_metrics_pipeline(
        model,
        d100k,
        output_path,
        delphi_labels,
        diseases_of_interest=None,
        filter_min_total=args.filter_min_total,
        disease_chunk_size=args.disease_chunk_size,
        device=device,
        seed=seed,
        n_bootstrap=args.n_bootstrap,
        threshold=args.threshold,
    )


if __name__ == "__main__":
    main()
