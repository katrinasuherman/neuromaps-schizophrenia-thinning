'''
env.py provides simple environment validation utilities for the Neuromaps
pipeline. It ensures that required external dependencies—particularly
the Connectome Workbench command-line tool (`wb_command`)—are available
and accessible in the system PATH.
'''

import os, subprocess

def ensure_workbench(path="/Applications/Workbench/bin_macosxub"):
    '''
    Verify that the Connectome Workbench (`wb_command`) executable is available
    and functioning. If found, returns its version string; otherwise raises an error.

    Notes
    - On macOS, Workbench is often installed under:
          /Applications/Workbench/bin_macosxub/
      On Linux, typical paths include:
          /usr/local/workbench/bin_linux64/
    - You can manually override this path when calling the function if your
      installation differs:
          >>> ensure_workbench("/usr/local/workbench/bin_linux64")
    - If Workbench is not strictly required (e.g., using pre-transformed maps),
      you can safely skip this check.
    '''
    os.environ["PATH"] = f"{path}:{os.environ.get('PATH','')}"
    try:
        cp = subprocess.run(["wb_command", "-version"], capture_output=True, text=True)
        cp.check_returncode()
        return cp.stdout.strip()
    except Exception as e:
        raise RuntimeError("Connectome Workbench not found on PATH") from e
