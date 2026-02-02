"""
DQN-based Spectrum Management for 5G-Radar Coexistence.

Implements:
  1. Standard DQN
  2. Double Dueling DQN with Prioritised Experience Replay (PER)

The agent observes a flattened state vector derived from the PRB graph and
selects a *joint* PRB allocation action (discretised).

State space  (from NRwithRADAR.m parameters):
  - CQI per UE per PRB  (4 × 25 = 100)
  - Radar active per PRB (25)
  - BLER per UE          (4)
  - Buffer status per UE  (4)
  Total = 133

Action space:
  Each PRB can be assigned to one of NUM_UES or left empty.
  Full combinatorial is too large → we use a *per-PRB* sequential approach
  or a reduced discrete action set.
"""

import random
import numpy as np
from collections import deque, namedtuple

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
NUM_RBS = 25
NUM_UES = 4
STATE_DIM = NUM_UES * NUM_RBS + NUM_RBS + NUM_UES + NUM_UES  # 133
ACTION_DIM = NUM_UES + 1  # per-PRB action: assign to UE 0-3 or leave empty

Transition = namedtuple("Transition", ("state", "action", "reward", "next_state", "done"))


# ---------------------------------------------------------------------------
# 1. Replay buffers
# ---------------------------------------------------------------------------

class ReplayBuffer:
    """Standard uniform replay buffer."""

    def __init__(self, capacity: int = 100_000):
        self.buffer = deque(maxlen=capacity)

    def push(self, *args):
        self.buffer.append(Transition(*args))

    def sample(self, batch_size: int):
        batch = random.sample(self.buffer, batch_size)
        return Transition(*zip(*batch))

    def __len__(self):
        return len(self.buffer)


class PrioritisedReplayBuffer:
    """Proportional prioritised experience replay (PER)."""

    def __init__(self, capacity: int = 100_000, alpha: float = 0.6):
        self.capacity = capacity
        self.alpha = alpha
        self.buffer = []
        self.priorities = np.zeros(capacity, dtype=np.float32)
        self.pos = 0

    def push(self, *args):
        max_prio = self.priorities[: len(self.buffer)].max() if self.buffer else 1.0
        if len(self.buffer) < self.capacity:
            self.buffer.append(Transition(*args))
        else:
            self.buffer[self.pos] = Transition(*args)
        self.priorities[self.pos] = max_prio
        self.pos = (self.pos + 1) % self.capacity

    def sample(self, batch_size: int, beta: float = 0.4):
        n = len(self.buffer)
        probs = self.priorities[:n] ** self.alpha
        probs /= probs.sum()
        indices = np.random.choice(n, batch_size, p=probs, replace=False)
        weights = (n * probs[indices]) ** (-beta)
        weights /= weights.max()
        batch = [self.buffer[i] for i in indices]
        return Transition(*zip(*batch)), indices, torch.tensor(weights, dtype=torch.float32)

    def update_priorities(self, indices, priorities):
        for idx, prio in zip(indices, priorities):
            self.priorities[idx] = abs(prio) + 1e-6

    def __len__(self):
        return len(self.buffer)


# ---------------------------------------------------------------------------
# 2. Network architectures
# ---------------------------------------------------------------------------

class DQNNetwork(nn.Module):
    """Standard DQN (multi-head: one Q-head per PRB)."""

    def __init__(self, state_dim: int = STATE_DIM, action_dim: int = ACTION_DIM,
                 hidden: int = 256):
        super().__init__()
        self.shared = nn.Sequential(
            nn.Linear(state_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
        )
        # One head per PRB → outputs (batch, NUM_RBS * ACTION_DIM)
        self.q_head = nn.Linear(hidden, NUM_RBS * action_dim)
        self.action_dim = action_dim

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        """Returns Q-values shaped (batch, NUM_RBS, ACTION_DIM)."""
        h = self.shared(state)
        q = self.q_head(h)
        return q.view(-1, NUM_RBS, self.action_dim)


class DuelingDQNNetwork(nn.Module):
    """Dueling architecture with separate value and advantage streams."""

    def __init__(self, state_dim: int = STATE_DIM, action_dim: int = ACTION_DIM,
                 hidden: int = 256):
        super().__init__()
        self.shared = nn.Sequential(
            nn.Linear(state_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
        )
        self.value_stream = nn.Sequential(
            nn.Linear(hidden, 128),
            nn.ReLU(),
            nn.Linear(128, NUM_RBS),  # V(s) per PRB
        )
        self.advantage_stream = nn.Sequential(
            nn.Linear(hidden, 128),
            nn.ReLU(),
            nn.Linear(128, NUM_RBS * action_dim),
        )
        self.action_dim = action_dim

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        h = self.shared(state)
        v = self.value_stream(h).unsqueeze(-1)  # (B, PRB, 1)
        a = self.advantage_stream(h).view(-1, NUM_RBS, self.action_dim)  # (B, PRB, A)
        q = v + a - a.mean(dim=-1, keepdim=True)
        return q


# ---------------------------------------------------------------------------
# 3. Spectrum environment (offline, data-driven)
# ---------------------------------------------------------------------------

class SpectrumEnv:
    """Gym-like environment that replays MATLAB simulation snapshots.

    Each step:
      - observation: CQI + radar + BLER + buffer status
      - action: per-PRB allocation vector
      - reward: weighted sum of throughput gain + radar avoidance + fairness
    """

    def __init__(self, dataset: list, reward_weights: dict | None = None):
        """
        Parameters
        ----------
        dataset : list of dicts with keys 'cqi', 'radar', 'bler', 'throughput', 'buffer'
        reward_weights : optional dict with 'throughput', 'radar_penalty', 'fairness'
        """
        self.dataset = dataset
        self.idx = 0
        self.rw = reward_weights or {
            "throughput": 1.0,
            "radar_penalty": -5.0,
            "fairness": 0.5,
        }

    def reset(self) -> np.ndarray:
        self.idx = 0
        return self._get_state()

    def step(self, action: np.ndarray):
        """
        action : (NUM_RBS,) int – UE assignment per PRB (0..NUM_UES, where NUM_UES = unallocated)
        """
        snap = self.dataset[self.idx]
        reward = self._compute_reward(action, snap)
        self.idx = (self.idx + 1) % len(self.dataset)
        done = self.idx == 0
        return self._get_state(), reward, done, {}

    def _get_state(self) -> np.ndarray:
        snap = self.dataset[self.idx % len(self.dataset)]
        cqi = snap["cqi"].flatten()  # (100,)
        radar = snap["radar"]  # (25,)
        bler = snap.get("bler", np.zeros(NUM_UES))
        buf = snap.get("buffer", np.zeros(NUM_UES))
        return np.concatenate([cqi, radar, bler, buf]).astype(np.float32)

    def _compute_reward(self, action: np.ndarray, snap: dict) -> float:
        cqi = snap["cqi"]  # (UE, PRB)
        radar = snap["radar"]  # (PRB,)

        # Throughput component: sum of CQI at allocated PRBs
        tp = 0.0
        for prb in range(NUM_RBS):
            ue = action[prb]
            if 0 <= ue < NUM_UES:
                tp += cqi[ue, prb]

        # Radar penalty: allocating on radar-active PRBs
        radar_pen = 0.0
        for prb in range(NUM_RBS):
            if radar[prb] > 0.5 and action[prb] < NUM_UES:
                radar_pen += 1.0

        # Fairness: Jain's index
        ue_alloc = np.zeros(NUM_UES)
        for prb in range(NUM_RBS):
            ue = action[prb]
            if 0 <= ue < NUM_UES:
                ue_alloc[ue] += cqi[ue, prb]
        sum_sq = ue_alloc.sum() ** 2
        sq_sum = (ue_alloc ** 2).sum()
        jain = sum_sq / (NUM_UES * sq_sum) if sq_sum > 0 else 0.0

        reward = (
            self.rw["throughput"] * tp / (NUM_RBS * 15)
            + self.rw["radar_penalty"] * radar_pen / NUM_RBS
            + self.rw["fairness"] * jain
        )
        return float(reward)


# ---------------------------------------------------------------------------
# 4. DQN Agent
# ---------------------------------------------------------------------------

class DQNAgent:
    """Double-Dueling DQN agent with optional PER."""

    def __init__(
        self,
        state_dim: int = STATE_DIM,
        action_dim: int = ACTION_DIM,
        hidden: int = 256,
        lr: float = 1e-4,
        gamma: float = 0.99,
        tau: float = 0.005,
        eps_start: float = 1.0,
        eps_end: float = 0.05,
        eps_decay: int = 10_000,
        use_per: bool = True,
        use_dueling: bool = True,
        buffer_size: int = 100_000,
        batch_size: int = 64,
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

        Net = DuelingDQNNetwork if use_dueling else DQNNetwork
        self.policy_net = Net(state_dim, action_dim, hidden).to(self.device)
        self.target_net = Net(state_dim, action_dim, hidden).to(self.device)
        self.target_net.load_state_dict(self.policy_net.state_dict())

        self.optimizer = optim.Adam(self.policy_net.parameters(), lr=lr)
        self.use_per = use_per
        if use_per:
            self.buffer = PrioritisedReplayBuffer(buffer_size)
        else:
            self.buffer = ReplayBuffer(buffer_size)

    @property
    def epsilon(self) -> float:
        return self.eps_end + (self.eps_start - self.eps_end) * np.exp(
            -self.steps / self.eps_decay
        )

    def select_action(self, state: np.ndarray) -> np.ndarray:
        """ε-greedy per-PRB action selection."""
        self.steps += 1
        if random.random() < self.epsilon:
            return np.random.randint(0, ACTION_DIM, size=NUM_RBS)
        with torch.no_grad():
            s = torch.tensor(state, dtype=torch.float32).unsqueeze(0).to(self.device)
            q = self.policy_net(s)  # (1, PRB, A)
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

        states = torch.tensor(np.array(batch.state), dtype=torch.float32, device=self.device)
        actions = torch.tensor(np.array(batch.action), dtype=torch.long, device=self.device)
        rewards = torch.tensor(np.array(batch.reward), dtype=torch.float32, device=self.device)
        next_states = torch.tensor(np.array(batch.next_state), dtype=torch.float32, device=self.device)
        dones = torch.tensor(np.array(batch.done), dtype=torch.float32, device=self.device)

        # Q(s,a) for each PRB
        q_vals = self.policy_net(states)  # (B, PRB, A)
        q_selected = q_vals.gather(2, actions.unsqueeze(-1)).squeeze(-1)  # (B, PRB)

        # Double DQN target
        with torch.no_grad():
            next_q_policy = self.policy_net(next_states)
            next_actions = next_q_policy.argmax(dim=2, keepdim=True)
            next_q_target = self.target_net(next_states)
            next_q = next_q_target.gather(2, next_actions).squeeze(-1)
            target = rewards.unsqueeze(1) + self.gamma * next_q * (1 - dones.unsqueeze(1))

        td_error = (q_selected - target).abs().mean(dim=1)
        loss = (weights * (q_selected - target).pow(2).mean(dim=1)).mean()

        self.optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(self.policy_net.parameters(), 10.0)
        self.optimizer.step()

        # Soft target update
        for tp, pp in zip(self.target_net.parameters(), self.policy_net.parameters()):
            tp.data.copy_(self.tau * pp.data + (1 - self.tau) * tp.data)

        if self.use_per:
            self.buffer.update_priorities(indices, td_error.detach().cpu().numpy())

        return loss.item()


# ---------------------------------------------------------------------------
# 5. Convenience: generate synthetic env data
# ---------------------------------------------------------------------------

def make_synthetic_env(num_snapshots: int = 5000, seed: int = 42) -> SpectrumEnv:
    rng = np.random.RandomState(seed)
    dataset = []
    for _ in range(num_snapshots):
        cqi = rng.randint(1, 16, size=(NUM_UES, NUM_RBS)).astype(np.float32)
        radar = np.zeros(NUM_RBS, dtype=np.float32)
        if rng.random() < 0.3:
            s = rng.randint(0, NUM_RBS)
            w = rng.randint(1, 8)
            radar[s: min(s + w, NUM_RBS)] = 1.0
        dataset.append({
            "cqi": cqi,
            "radar": radar,
            "bler": rng.uniform(0, 0.3, NUM_UES).astype(np.float32),
            "buffer": rng.uniform(0, 1, NUM_UES).astype(np.float32),
        })
    return SpectrumEnv(dataset)
