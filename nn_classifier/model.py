"""
MLP for P(homolog | features).
"""
from __future__ import annotations
import torch
import torch.nn as nn


class HomologyClassifier(nn.Module):
    def __init__(self, in_dim: int, hidden: tuple = (128, 128), dropout: float = 0.1):
        super().__init__()
        layers = []
        prev = in_dim
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # Returns logits — apply sigmoid downstream for probability.
        return self.net(x).squeeze(-1)
