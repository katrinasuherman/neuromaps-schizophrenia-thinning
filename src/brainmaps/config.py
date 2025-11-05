# src/brainmaps/config.py

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


'''
config.py defines the Config class, which holds project-wide settings
such as output directories, random seeds, and the number of permutations
for spin tests. It also provides a method to load these parameters
from a JSON configuration file.
'''

@dataclass
class Config:
    out_dir: Path = Path("out")
    seed: int = 42
    n_perm: int = 1000

    @staticmethod
    def from_json(path: Optional[str] = "configs/config.json"):
        '''
        Load configuration values from a JSON file and return a Config object.

        :param path: Path to the JSON configuration file. If None, default values are used.
        :return: Config object with attributes populated from the JSON file (if found),
                 or defaults otherwise.
        '''
        if not path or not Path(path).exists():
            print(f"[config] using default settings (no config file found at {path})")
            return Config()

        with open(path, "r") as f:
            data = json.load(f)

        out_dir = Path(data.get("out_dir", "out"))
        seed = int(data.get("seed", 42))
        n_perm = int(data.get("n_perm", 1000))

        print(f"[config] loaded: out_dir={out_dir}, seed={seed}, n_perm={n_perm}")
        return Config(out_dir=out_dir, seed=seed, n_perm=n_perm)


def save_default_config(path: str = "configs/config.json"):
    '''
    Create a default configuration file (if missing) with sensible defaults.

    :param path: Path where the default config file should be saved.
    :return: None
    '''
    cfg_path = Path(path)
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    default_data = {
        "out_dir": "out",
        "seed": 42,
        "n_perm": 1000
    }
    with open(cfg_path, "w") as f:
        json.dump(default_data, f, indent=2)
    print(f"[config] default config saved to {cfg_path}")
