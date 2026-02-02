"""
Hybrid GNN-DQN for 5G-Radar Spectrum Management.

Architecture:
    PRB Graph  ──►  GNN Feature Extractor  ──►  graph embedding (128-d)
                                                       │
                                        ┌──────────────┘
                                        ▼
                    [graph_emb ‖ flat_state]  ──►  Dueling DQN  ──►  per-PRB Q-values

The GNN provides structural awareness of PRB relationships and radar
interference patterns, while the DQN handles sequential decision-making
with exploration.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from ml.gnn_spectrum import GNNFeatureExtractor, NODE_FEATURE_DIM
from ml.dqn_spectrum import (
    NUM_RBS,
    NUM_UES,
    ACTION_DIM,
    STATE_DIM,
    ReplayBuffer,
    PrioritisedReplayBuffer,
    SpectrumEnv,
    Transition,
)
from ml.data_pipeline import create_prb_graph, compute_radar_overlap

# ---------------------------------------------------------------------------
# Hybrid network
# ---------------------------------------------------------------------------

class HybridGNNDQN(nn.Module):
    """Combines GNN graph-level embedding with Dueling DQN."""

    def __init__(
        self,
        gnn_out_dim: int = 128,
        state_dim: int = STATE_DIM,
        action_dim: int = ACTION_DIM,
        hidden: int = 256,
    ):
        super().__init__()
        self.gnn = GNNFeatureExtractor(
            in_channels=NODE_FEATURE_DIM,
            hidden=64,
            out_dim=gnn_out_dim,
            num_layers=3,
        )
        combined_dim = gnn_out_dim + state_dim

        self.shared = nn.Sequential(
            nn.Linear(combined_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
        )
        self.value = nn.Sequential(
            nn.Linear(hidden, 128),
            nn.ReLU(),
            nn.Linear(128, NUM_RBS),
        )
        self.advantage = nn.Sequential(
            nn.Linear(hidden, 128),
            nn.ReLU(),
            nn.Linear(128, NUM_RBS * action_dim),
        )
        self.action_dim = action_dim

    def forward(self, graph_data, flat_state: torch.Tensor) -> torch.Tensor:
        """
        Parameters
        ----------
        graph_data : PyG Data (batched)
        flat_state : (batch, STATE_DIM)

        Returns
        -------
        Q-values : (batch, NUM_RBS, ACTION_DIM)
        """
        graph_emb = self.gnn(graph_data)  # (batch, gnn_out_dim)
        combined = torch.cat([graph_emb, flat_state], dim=1)
        h = self.shared(combined)
        v = self.value(h).unsqueeze(-1)
        a = self.advantage(h).view(-1, NUM_RBS, self.action_dim)
        return v + a - a.mean(dim=-1, keepdim=True)


# ---------------------------------------------------------------------------
# Hybrid agent
# ---------------------------------------------------------------------------

class HybridAgent:
    """Hybrid GNN-DQN agent combining graph perception with RL."""

    def __init__(
        self,
        lr: float = 3e-4,
        gamma: float = 0.99,
        tau: float = 0.005,
        eps_start: float = 1.0,
        eps_end: float = 0.05,
        eps_decay: int = 10_000,
        use_per: bool = True,
        buffer_size: int = 100_000,
        batch_size: int = 32,
        device: str = "cpu",
    ):
        self.device = torch.device(device)
        self.gamma = gamma
        self.tau = tau
        self.eps_start = eps_start
        self.eps_end = eps_end
        self.eps_decay = eps_decay
        self.batch_size = batch_size
        self.steps = 0

        self.policy_net = HybridGNNDQN().to(self.device)
        self.target_net = HybridGNNDQN().to(self.device)
        self.target_net.load_state_dict(self.policy_net.state_dict())

        self.optimizer = optim.Adam(self.policy_net.parameters(), lr=lr)
        self.use_per = use_per
        self.buffer = PrioritisedReplayBuffer(buffer_size) if use_per else ReplayBuffer(buffer_size)

    @property
    def epsilon(self) -> float:
        return self.eps_end + (self.eps_start - self.eps_end) * np.exp(
            -self.steps / self.eps_decay
        )

    def _state_to_graph(self, state: np.ndarray):
        """Convert flat state vector back to a PyG graph for the GNN."""
        cqi = state[: NUM_UES * NUM_RBS].reshape(NUM_UES, NUM_RBS)
        radar = state[NUM_UES * NUM_RBS: NUM_UES * NUM_RBS + NUM_RBS]
        bler = state[NUM_UES * NUM_RBS + NUM_RBS: NUM_UES * NUM_RBS + NUM_RBS + NUM_UES]
        graph = create_prb_graph(cqi, radar, bler)
        graph.batch = torch.zeros(graph.x.size(0), dtype=torch.long)
        return graph

    def select_action(self, state: np.ndarray) -> np.ndarray:
        self.steps += 1
        if np.random.random() < self.epsilon:
            return np.random.randint(0, ACTION_DIM, size=NUM_RBS)

        with torch.no_grad():
            graph = self._state_to_graph(state).to(self.device)
            flat = torch.tensor(state, dtype=torch.float32).unsqueeze(0).to(self.device)
            q = self.policy_net(graph, flat)
            return q.squeeze(0).argmax(dim=1).cpu().numpy()

    def store(self, state, action, reward, next_state, done):
        self.buffer.push(state, action, reward, next_state, done)

    def update(self) -> float:
        if len(self.buffer) < self.batch_size:
            return 0.0

        if self.use_per:
            batch, indices, weights = self.buffer.sample(self.batch_size)
            weights = weights.to(self.device)
        else:
            batch = self.buffer.sample(self.batch_size)
            weights = torch.ones(self.batch_size, device=self.device)

        states = np.array(batch.state)
        actions = torch.tensor(np.array(batch.action), dtype=torch.long, device=self.device)
        rewards = torch.tensor(np.array(batch.reward), dtype=torch.float32, device=self.device)
        next_states = np.array(batch.next_state)
        dones = torch.tensor(np.array(batch.done), dtype=torch.float32, device=self.device)

        # Build batched graphs
        from torch_geometric.data import Batch
        graphs = Batch.from_data_list([self._state_to_graph(s) for s in states]).to(self.device)
        next_graphs = Batch.from_data_list([self._state_to_graph(s) for s in next_states]).to(self.device)

        flat_s = torch.tensor(states, dtype=torch.float32, device=self.device)
        flat_ns = torch.tensor(next_states, dtype=torch.float32, device=self.device)

        q_vals = self.policy_net(graphs, flat_s)
        q_selected = q_vals.gather(2, actions.unsqueeze(-1)).squeeze(-1)

        with torch.no_grad():
            next_q_policy = self.policy_net(next_graphs, flat_ns)
            next_actions = next_q_policy.argmax(dim=2, keepdim=True)
            next_q_target = self.target_net(next_graphs, flat_ns)
            next_q = next_q_target.gather(2, next_actions).squeeze(-1)
            target = rewards.unsqueeze(1) + self.gamma * next_q * (1 - dones.unsqueeze(1))

        td_error = (q_selected - target).abs().mean(dim=1)
        loss = (weights * (q_selected - target).pow(2).mean(dim=1)).mean()

        self.optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(self.policy_net.parameters(), 10.0)
        self.optimizer.step()

        for tp, pp in zip(self.target_net.parameters(), self.policy_net.parameters()):
            tp.data.copy_(self.tau * pp.data + (1 - self.tau) * tp.data)

        if self.use_per:
            self.buffer.update_priorities(indices, td_error.detach().cpu().numpy())

        return loss.item()

    def save(self, path: str):
        torch.save({
            "policy_net": self.policy_net.state_dict(),
            "target_net": self.target_net.state_dict(),
            "optimizer": self.optimizer.state_dict(),
            "steps": self.steps,
        }, path)

    def load(self, path: str):
        ckpt = torch.load(path, map_location=self.device)
        self.policy_net.load_state_dict(ckpt["policy_net"])
        self.target_net.load_state_dict(ckpt["target_net"])
        self.optimizer.load_state_dict(ckpt["optimizer"])
        self.steps = ckpt["steps"]
