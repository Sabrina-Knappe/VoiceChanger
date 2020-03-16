#this code is to test the pysndfx package
infile = 'my_audio_file.wav'
outfile = 'my_processed_audio_file.ogg'
import pysndfx

fx = (
    AudioEffectsChain()
    .highshelf()
    .reverb()
    .phaser()
    .delay()
    .lowshelf()
)
# Apply phaser and reverb directly to an audio file.
fx(infile, outfile)

# Or, apply the effects directly to a ndarray.
from librosa import load
y, sr = load(infile, sr=None)
y = fx(y)

# Apply the effects and return the results as a ndarray.
y = fx(infile)

# Apply the effects to a ndarray but store the resulting audio to disk.
fx(x, outfile)