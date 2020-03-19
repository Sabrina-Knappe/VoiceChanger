import pysndfx
from pysndfx import AudioEffectsChain
import librosa
from librosa import load
import pyo

#this code is to test the pysndfx package
infile = 'audio/273175__xserra__la-vaca-cega-xavier.wav'
outfile = 'audio/my_processed_audio_file.wav'


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
