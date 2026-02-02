"""
GNN-based Spectrum Management for 5G-Radar Coexistence.

Three architectures:
  1. GCN  – Graph Convolutional Network
  2. GAT  – Graph Attention Network
  3. SAGE – GraphSAGE

Each model maps a PRB graph (from data_pipeline) to per-PRB allocation
decisions (NUM_UES + 1 classes: one per UE + "do-not-allocate").

Designed to match the simulation parameters in NRwithRADAR.m:
  NUM_RBS = 25, NUM_UES = 4, node features = 10
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    from torch_geometric.nn import GCNConv, GATConv, SAGEConv, global_mean_pool
    from torch_geometric.data import Data, DataLoader as PyGDataLoader
except ImportError:
    raise ImportError(
        "torch_geometric is required. Install via:\n"
        "  pip install torch-geometric torch-scatter torch-sparse"
    )

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
NUM_RBS = 25
NUM_UES = 4
NODE_FEATURE_DIM = 10  # from data_pipeline.create_prb_graph
NUM_CLASSES = NUM_UES + 1  # 4 UEs + 1 "no allocate"


# ---------------------------------------------------------------------------
# 1. GCN Model
# ---------------------------------------------------------------------------
class GCNSpectrum(nn.Module):
    """Graph Convolutional Network for PRB allocation."""

    def __init__(
        self,
        in_channels: int = NODE_FEATURE_DIM,
        hidden: int = 64,
        num_layers: int = 3,
        out_channels: int = NUM_CLASSES,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        self.convs.append(GCNConv(in_channels, hidden))
        self.bns.append(nn.BatchNorm1d(hidden))
        for _ in range(num_layers - 2):
            self.convs.append(GCNConv(hidden, hidden))
            self.bns.append(nn.BatchNorm1d(hidden))
        self.convs.append(GCNConv(hidden, hidden))
        self.bns.append(nn.BatchNorm1d(hidden))
        self.head = nn.Linear(hidden, out_channels)
        self.dropout = dropout

    def forward(self, data: Data) -> torch.Tensor:
        x, edge_index = data.x, data.edge_index
        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        return self.head(x)  # (num_nodes, NUM_CLASSES)


# ---------------------------------------------------------------------------
# 2. GAT Model
# ---------------------------------------------------------------------------
class GATSpectrum(nn.Module):
    """Graph Attention Network for PRB allocation."""

    def __init__(
        self,
        in_channels: int = NODE_FEATURE_DIM,
        hidden: int = 64,
        num_layers: int = 3,
        out_channels: int = NUM_CLASSES,
        heads: int = 4,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        self.convs.append(GATConv(in_channels, hidden, heads=heads, concat=False))
        self.bns.append(nn.BatchNorm1d(hidden))
        for _ in range(num_layers - 2):
            self.convs.append(GATConv(hidden, hidden, heads=heads, concat=False))
            self.bns.append(nn.BatchNorm1d(hidden))
        self.convs.append(GATConv(hidden, hidden, heads=heads, concat=False))
        self.bns.append(nn.BatchNorm1d(hidden))
        self.head = nn.Linear(hidden, out_channels)
        self.dropout = dropout

    def forward(self, data: Data) -> torch.Tensor:
        x, edge_index = data.x, data.edge_index
        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        return self.head(x)


# ---------------------------------------------------------------------------
# 3. GraphSAGE Model
# ---------------------------------------------------------------------------
class SAGESpectrum(nn.Module):
    """GraphSAGE for PRB allocation."""

    def __init__(
        self,
        in_channels: int = NODE_FEATURE_DIM,
        hidden: int = 64,
        num_layers: int = 3,
        out_channels: int = NUM_CLASSES,
        dropout: float = 0.2,
    ):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        self.convs.append(SAGEConv(in_channels, hidden))
        self.bns.append(nn.BatchNorm1d(hidden))
        for _ in range(num_layers - 2):
            self.convs.append(SAGEConv(hidden, hidden))
            self.bns.append(nn.BatchNorm1d(hidden))
        self.convs.append(SAGEConv(hidden, hidden))
        self.bns.append(nn.BatchNorm1d(hidden))
        self.head = nn.Linear(hidden, out_channels)
        self.dropout = dropout

    def forward(self, data: Data) -> torch.Tensor:
        x, edge_index = data.x, data.edge_index
        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        return self.head(x)


# ---------------------------------------------------------------------------
# 4. Graph-level feature extractor (used by Hybrid GNN-DQN)
# ---------------------------------------------------------------------------
class GNNFeatureExtractor(nn.Module):
    """Extracts a fixed-size graph-level embedding from a PRB graph.

    Used as the perception backbone in the Hybrid GNN-DQN agent.
    """

    def __init__(
        self,
        in_channels: int = NODE_FEATURE_DIM,
        hidden: int = 64,
        out_dim: int = 128,
        num_layers: int = 3,
    ):
        super().__init__()
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        ch = in_channels
        for _ in range(num_layers):
            self.convs.append(GCNConv(ch, hidden))
            self.bns.append(nn.BatchNorm1d(hidden))
            ch = hidden
        self.fc = nn.Linear(hidden, out_dim)

    def forward(self, data: Data) -> torch.Tensor:
        x, edge_index, batch = data.x, data.edge_index, data.batch
        for conv, bn in zip(self.convs, self.bns):
            x = F.relu(bn(conv(x, edge_index)))
        # Global mean pooling → graph-level
        x = global_mean_pool(x, batch)
        return self.fc(x)  # (batch_size, out_dim)


# ---------------------------------------------------------------------------
# 5. Training utilities
# ---------------------------------------------------------------------------

def train_epoch(model, loader, optimizer, device="cpu"):
    model.train()
    total_loss = 0.0
    correct = 0
    total = 0
    criterion = nn.CrossEntropyLoss(ignore_index=-1)

    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        out = model(batch)  # (N, C)
        loss = criterion(out, batch.y)
        loss.backward()
        optimizer.step()

        total_loss += loss.item() * batch.num_nodes
        pred = out.argmax(dim=1)
        mask = batch.y >= 0
        correct += (pred[mask] == batch.y[mask]).sum().item()
        total += mask.sum().item()

    return total_loss / max(total, 1), correct / max(total, 1)


@torch.no_grad()
def eval_epoch(model, loader, device="cpu"):
    model.eval()
    total_loss = 0.0
    correct = 0
    total = 0
    criterion = nn.CrossEntropyLoss(ignore_index=-1)

    for batch in loader:
        batch = batch.to(device)
        out = model(batch)
        loss = criterion(out, batch.y)
        total_loss += loss.item() * batch.num_nodes
        pred = out.argmax(dim=1)
        mask = batch.y >= 0
        correct += (pred[mask] == batch.y[mask]).sum().item()
        total += mask.sum().item()

    return total_loss / max(total, 1), correct / max(total, 1)


def get_model(name: str, **kwargs) -> nn.Module:
    """Factory: return a GNN model by name."""
    models = {
        "gcn": GCNSpectrum,
        "gat": GATSpectrum,
        "sage": SAGESpectrum,
    }
    if name.lower() not in models:
        raise ValueError(f"Unknown model: {name}. Choose from {list(models)}")
    return models[name.lower()](**kwargs)
