import os, subprocess

def ensure_workbench(path="/Applications/Workbench/bin_macosxub"):
    os.environ["PATH"] = f"{path}:{os.environ.get('PATH','')}"
    try:
        cp = subprocess.run(["wb_command", "-version"], capture_output=True, text=True)
        cp.check_returncode()
        return cp.stdout.strip()
    except Exception as e:
        raise RuntimeError("Connectome Workbench not found on PATH") from e
