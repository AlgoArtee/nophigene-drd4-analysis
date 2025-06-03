# ---- 1a. Base Image: Python 3.10-slim (Debian 12 “bullseye”) ----
FROM python:3.10-slim

# Explanation:
# 1. We choose python:3.10-slim because:
#    - DeepChem requires Python < 3.11 (so 3.10 is the latest compatible). :contentReference[oaicite:9]{index=9}
#    - Matplotlib and other scientific libraries are fully supported on 3.10. 
#    - Slim variants minimize image size by excluding extra Debian components. :contentReference[oaicite:11]{index=11}

# ---- 1b. Patch All System Packages to Mitigate High-Severity CVEs ----
RUN apt-get update && \
    apt-get upgrade -y && \
    rm -rf /var/lib/apt/lists/*

# Explanation:
# 1. Debian 12 “bullseye” (used by python:3.10-slim) originally contained:
#    - CVE-2024-45491 (libexpat heap overflow, CVSS 9.8)
#    - CVE-2024-45492 (libexpat stack buffer overflow, CVSS 9.8)
#    - CVE-2024-5535  (glibc information exposure, CVSS 9.1)
#    :contentReference[oaicite:12]{index=12}
# 2. Running 'apt-get upgrade' pulls in the backported fixes (e.g., libexpat3 2.4.9-1+deb12uX, glibc 2.31-13+deb11uX), eliminating high CVEs. :contentReference[oaicite:13]{index=13}
# 3. Cleaning '/var/lib/apt/lists/' reduces image size and prevents stale indices. :contentReference[oaicite:14]{index=14}

# ---- 1c. Install Build-Time Dependencies Needed for Bioinformatics C Extensions ----
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    wget \
    # (Add libpq-dev, git, wget, etc., if needed later)
    && rm -rf /var/lib/apt/lists/*

# Explanation:
# 1. 'build-essential' provides gcc, g++, make, etc., required to compile C extensions. :contentReference[oaicite:15]{index=15}
# 2. 'libbz2-dev' and 'liblzma-dev' enable building modules that rely on bzip2/xz libraries (e.g., pysam, scikit-allel). :contentReference[oaicite:16]{index=16}
# 3. 'libcurl4-openssl-dev' and 'libssl-dev' are needed for TLS and HTTP operations in packages like RDKit or PyCBC. :contentReference[oaicite:17]{index=17}
# 4. Using '--no-install-recommends' avoids extra packages, slimming the image. :contentReference[oaicite:18]{index=18}

# ---- 1d. Create a Non-Root User for Better Security (Optional but Best Practice) ----
ARG USERNAME=appuser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME && \
    useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

# Explanation:
# 1. Creating a non-root user avoids running Python scripts as root inside the container—a security best practice. :contentReference[oaicite:19]{index=19}
# 2. Debian ’bullseye’ will add this user without extra system packages. :contentReference[oaicite:20]{index=20}

# ---- 1e. Switch to the Non-Root User and Set Working Directory ----
USER $USERNAME
WORKDIR /home/$USERNAME/app

# Explanation:
# 1. All subsequent 'COPY' and 'RUN' steps happen as 'appuser,' minimizing root-owned files. :contentReference[oaicite:21]{index=21}
# 2. The working dir '/home/appuser/app' becomes the project root inside the container. :contentReference[oaicite:22]{index=22}

# 2. Copy the requirements file first to leverage Docker layer caching ----
COPY --chown=$USERNAME:$USERNAME requirements.txt .

# ---- 2a. Install Python dependencies inside the container ----
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Explanation:
# - Inside a Linux container, 'pip install pysam==0.23.1' fetches a prebuilt wheel, avoiding compile errors. :contentReference[oaicite:13]{index=13}

# ---- 2b. Copy the rest of the source code into the image ----
COPY --chown=$USERNAME:$USERNAME . .

# ---- 2c. Define the default command to run your analysis script ----
#     Modify arguments as needed once your CLI options are finalized.
ENTRYPOINT ["python", "src/analysis.py"]

# Explanation:
# - When the container runs, it executes 'python src/analysis.py'. :contentReference[oaicite:14]{index=14}