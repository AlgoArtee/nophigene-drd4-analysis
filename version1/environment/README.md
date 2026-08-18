# Retired Version 1 Python environment

The former root `.venv` occupied approximately 575 MB and referenced a removed Python 3.10 installation at `C:\Users\Mewxy\AppData\Local\Programs\Python\Python310\python.exe`. Its interpreter, `pip`, and dependency checks could no longer start, so the environment was non-portable and unusable.

The directory was not archived byte-for-byte because Python virtual environments are generated caches containing absolute interpreter paths. Recreate a historical environment only if needed by installing the archived `requirements-research.txt` together with the appropriate historical application requirements in a fresh Python 3.10 environment.

Version 2 uses the root `.venv-v2` development/test environment (Python 3.12) and the pinned Linux container for supported execution.
