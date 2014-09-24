from matplotlib.colors import ListedColormap
from .data import Spectrogram, PCM, gaussian_spectrum, strfpak_spectrum
from .file import Wave, read, si, si_from_name

# This munge necessary to deal with scikits.samplerate
# This makes "namespace packages" available.
__import__('pkg_resources').declare_namespace(__name__)