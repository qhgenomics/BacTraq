import warnings

# Silences pandas FutureWarnings triggered by groupby/apply patterns used throughout
# clustering.py and cluster_trace.py. Set at package import time so it's in effect
# regardless of which submodule is imported first.
warnings.simplefilter(action='ignore', category=FutureWarning)

__version__ = "1.1.0-rc.2"
