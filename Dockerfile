FROM --platform=linux/amd64 python:3.11

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# ---------------------------------------------------------------------
# System dependencies for neuroimaging & plotting
# ---------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender1 \
    libglu1-mesa \
    libgl1-mesa-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------
# Install Connectome Workbench
# ---------------------------------------------------------------------
COPY workbench /opt/workbench

# Fix permissions on binaries (otherwise Docker can't execute them)
RUN chmod -R a+x /opt/workbench/bin_linux64 \
    && chmod -R a+x /opt/workbench/exe_linux64

# Add Workbench bin folder to PATH
ENV PATH="/opt/workbench/bin_linux64:${PATH}"

# ---------------------------------------------------------------------
# Install Python dependencies
# ---------------------------------------------------------------------
WORKDIR /app

# Copy requirements first (for layer caching)
COPY requirements.txt /app/

RUN pip install --upgrade pip setuptools wheel && \
    pip install -r requirements.txt

# ---------------------------------------------------------------------
# Copy project code
# ---------------------------------------------------------------------
COPY . /app

# Default interactive shell
CMD ["bash"]
