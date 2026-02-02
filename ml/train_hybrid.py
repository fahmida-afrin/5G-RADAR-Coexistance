"""
Training script for Hybrid GNN-DQN and standalone DQN agents.

Usage:
    python -m ml.train_hybrid --agent hybrid --episodes 1000
    python -m ml.train_hybrid --agent dqn --episodes 1000
"""

import os
import argparse
import json
import time
import numpy as np

import torch

from ml.dqn_spectrum import DQNAgent, make_synthetic_env, NUM_RBS
from ml.hybrid_gnn_dqn import HybridAgent


def train(args):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    # -- Environment --
    env = make_synthetic_env(num_snapshots=args.env_size, seed=args.seed)
    eval_env = make_synthetic_env(num_snapshots=500, seed=args.seed + 1)

    # -- Agent --
    if args.agent == "dqn":
        agent = DQNAgent(
            hidden=args.hidden,
            lr=args.lr,
            gamma=args.gamma,
            eps_decay=args.eps_decay,
            use_per=args.per,
            use_dueling=args.dueling,
            batch_size=args.batch_size,
            device=device,
        )
    elif args.agent == "hybrid":
        agent = HybridAgent(
            lr=args.lr,
            gamma=args.gamma,
            eps_decay=args.eps_decay,
            use_per=args.per,
            batch_size=args.batch_size,
            device=device,
        )
    else:
        raise ValueError(f"Unknown agent: {args.agent}")

    print(f"Agent: {args.agent}")

    os.makedirs(args.save_dir, exist_ok=True)
    history = {"episode_reward": [], "episode_loss": [], "epsilon": [], "eval_reward": []}
    best_eval = -np.inf

    for episode in range(1, args.episodes + 1):
        state = env.reset()
        total_reward = 0.0
        total_loss = 0.0
        steps = 0
        done = False

        while not done:
            action = agent.select_action(state)
            next_state, reward, done, _ = env.step(action)
            agent.store(state, action, reward, next_state, done)
            loss = agent.update()

            total_reward += reward
            total_loss += loss
            steps += 1
            state = next_state

            if steps > args.max_steps:
                break

        history["episode_reward"].append(total_reward)
        history["episode_loss"].append(total_loss / max(steps, 1))
        history["epsilon"].append(agent.epsilon)

        # Evaluation
        if episode % args.eval_every == 0:
            eval_reward = _evaluate_agent(agent, eval_env, num_episodes=5)
            history["eval_reward"].append(eval_reward)

            if eval_reward > best_eval:
                best_eval = eval_reward
                if hasattr(agent, "save"):
                    agent.save(os.path.join(args.save_dir, f"best_{args.agent}.pt"))
                else:
                    torch.save(
                        agent.policy_net.state_dict(),
                        os.path.join(args.save_dir, f"best_{args.agent}.pt"),
                    )

            if episode % args.log_every == 0:
                print(
                    f"Ep {episode:5d}/{args.episodes} | "
                    f"Reward: {total_reward:7.2f} | "
                    f"Eval: {eval_reward:7.2f} | "
                    f"Loss: {total_loss / max(steps, 1):.4f} | "
                    f"ε: {agent.epsilon:.3f}"
                )

    # Save final
    if hasattr(agent, "save"):
        agent.save(os.path.join(args.save_dir, f"final_{args.agent}.pt"))
    else:
        torch.save(
            agent.policy_net.state_dict(),
            os.path.join(args.save_dir, f"final_{args.agent}.pt"),
        )

    with open(os.path.join(args.save_dir, f"history_{args.agent}.json"), "w") as f:
        json.dump(history, f, indent=2)

    print(f"\nTraining complete. Best eval reward: {best_eval:.2f}")
    print(f"Checkpoints saved to {args.save_dir}")


def _evaluate_agent(agent, env, num_episodes: int = 5) -> float:
    """Run agent greedily on eval env."""
    total = 0.0
    orig_steps = agent.steps
    for _ in range(num_episodes):
        state = env.reset()
        done = False
        steps = 0
        while not done and steps < 500:
            # Force greedy
            old_eps = agent.eps_start
            agent.eps_start = 0.0
            agent.eps_end = 0.0
            action = agent.select_action(state)
            agent.eps_start = old_eps
            agent.eps_end = 0.05

            state, reward, done, _ = env.step(action)
            total += reward
            steps += 1
    agent.steps = orig_steps  # restore
    return total / num_episodes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train DQN/Hybrid agent")
    parser.add_argument("--agent", type=str, default="hybrid", choices=["dqn", "hybrid"])
    parser.add_argument("--episodes", type=int, default=1000)
    parser.add_argument("--max-steps", type=int, default=500)
    parser.add_argument("--lr", type=float, default=3e-4)
    parser.add_argument("--gamma", type=float, default=0.99)
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--eps-decay", type=int, default=10000)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--per", action="store_true", default=True)
    parser.add_argument("--dueling", action="store_true", default=True)
    parser.add_argument("--env-size", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--save-dir", type=str, default="data/results/rl_checkpoints")
    parser.add_argument("--log-every", type=int, default=10)
    parser.add_argument("--eval-every", type=int, default=10)
    args = parser.parse_args()

    train(args)
