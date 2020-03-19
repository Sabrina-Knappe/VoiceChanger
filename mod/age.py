import pysndfx
from pysndfx import AudioEffectsChain
import librosa
from librosa import load
import numpy as np 


def age(infile, outfile, value):
    fx = (
        AudioEffectsChain()
        .pitch(shift=10)
        # .tremolo(freq=3000)
        # .speed(factor=1.1)
    )
    # Apply phaser and reverb directly to an audio file.
    fx(infile, outfile)
    return outfile
