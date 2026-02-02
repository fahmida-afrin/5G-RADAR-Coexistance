"""
Training script for GNN-based PRB allocation models.

Usage:
    python -m ml.train_gnn --model gcn --epochs 200 --lr 1e-3
    python -m ml.train_gnn --model gat --epochs 200
    python -m ml.train_gnn --model sage --epochs 200
"""

import os
import argparse
import json
import time
import numpy as np

import torch
import torch.optim as optim
from torch_geometric.loader import DataLoader as PyGDataLoader

from ml.gnn_spectrum import get_model, train_epoch, eval_epoch
from ml.data_pipeline import generate_synthetic_dataset, split_dataset


def train(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # -- Data --
    data_path = args.data
    if data_path and os.path.exists(data_path):
        print(f"Loading dataset from {data_path}")
        dataset = torch.load(data_path)
    else:
        print(f"Generating synthetic dataset ({args.samples} samples)...")
        dataset = generate_synthetic_dataset(
            num_samples=args.samples,
            radar_duty_cycle=args.radar_dc,
        )

    train_data, test_data = split_dataset(dataset, train_ratio=0.8)
    print(f"Train: {len(train_data)}, Test: {len(test_data)}")

    train_loader = PyGDataLoader(train_data, batch_size=args.batch_size, shuffle=True)
    test_loader = PyGDataLoader(test_data, batch_size=args.batch_size, shuffle=False)

    # -- Model --
    model = get_model(
        args.model,
        hidden=args.hidden,
        num_layers=args.layers,
        dropout=args.dropout,
    ).to(device)

    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    print(f"Model: {args.model.upper()}, params: {sum(p.numel() for p in model.parameters()):,}")

    # -- Training loop --
    os.makedirs(args.save_dir, exist_ok=True)
    history = {"train_loss": [], "train_acc": [], "test_loss": [], "test_acc": []}
    best_acc = 0.0

    for epoch in range(1, args.epochs + 1):
        t0 = time.time()
        train_loss, train_acc = train_epoch(model, train_loader, optimizer, device)
        test_loss, test_acc = eval_epoch(model, test_loader, device)
        scheduler.step()
        elapsed = time.time() - t0

        history["train_loss"].append(train_loss)
        history["train_acc"].append(train_acc)
        history["test_loss"].append(test_loss)
        history["test_acc"].append(test_acc)

        if epoch % args.log_every == 0:
            print(
                f"Epoch {epoch:4d}/{args.epochs} | "
                f"Train Loss: {train_loss:.4f} Acc: {train_acc:.4f} | "
                f"Test Loss: {test_loss:.4f} Acc: {test_acc:.4f} | "
                f"{elapsed:.1f}s"
            )

        if test_acc > best_acc:
            best_acc = test_acc
            torch.save(model.state_dict(), os.path.join(args.save_dir, f"best_{args.model}.pt"))

        if epoch % args.ckpt_every == 0:
            torch.save(
                model.state_dict(),
                os.path.join(args.save_dir, f"{args.model}_epoch{epoch}.pt"),
            )

    # Save final
    torch.save(model.state_dict(), os.path.join(args.save_dir, f"final_{args.model}.pt"))
    with open(os.path.join(args.save_dir, f"history_{args.model}.json"), "w") as f:
        json.dump(history, f, indent=2)

    print(f"\nTraining complete. Best test accuracy: {best_acc:.4f}")
    print(f"Checkpoints saved to {args.save_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train GNN spectrum model")
    parser.add_argument("--model", type=str, default="gcn", choices=["gcn", "gat", "sage"])
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--wd", type=float, default=1e-5)
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--layers", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--samples", type=int, default=2000)
    parser.add_argument("--radar-dc", type=float, default=0.07)
    parser.add_argument("--data", type=str, default="data/processed/dataset.pt")
    parser.add_argument("--save-dir", type=str, default="data/results/gnn_checkpoints")
    parser.add_argument("--log-every", type=int, default=10)
    parser.add_argument("--ckpt-every", type=int, default=50)
    args = parser.parse_args()

    train(args)
