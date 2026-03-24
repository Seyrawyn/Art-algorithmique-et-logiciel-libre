from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field


@dataclass(slots=True)
class LowPass:
    alpha: float
    value: float = 0.0

    def step(self, sample: float) -> float:
        self.value += (sample - self.value) * self.alpha
        return self.value


@dataclass(slots=True)
class Envelope:
    attack: float
    release: float
    value: float = 0.0

    def step(self, sample: float) -> float:
        rate = self.attack if sample > self.value else self.release
        self.value += (sample - self.value) * rate
        return self.value


@dataclass(slots=True)
class SignalChain:
    taps: deque[float] = field(default_factory=lambda: deque(maxlen=12))
    low: LowPass = field(default_factory=lambda: LowPass(alpha=0.14))
    env: Envelope = field(default_factory=lambda: Envelope(attack=0.35, release=0.08))

    def push(self, sample: float) -> float:
        self.taps.append(sample)
        average = sum(self.taps) / len(self.taps)
        filtered = self.low.step(average)
        return self.env.step(filtered)


def pulse(index: int, t: float) -> float:
    phase = (index * 0.23) + (t * 2.1)
    return 0.5 + 0.5 * math.sin(phase + math.sin(phase * 0.7))


if __name__ == "__main__":
    chain = SignalChain()
    output = []
    for frame in range(48):
        t = frame / 24.0
        sample = pulse(frame, t)
        output.append(round(chain.push(sample), 4))
    print(output)
