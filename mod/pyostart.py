import pyo

snds_path = "audio"
s = pyo.Server().boot()
s.start()
a = pyo.SfPlayer(snds_path+"/my_processed_audio_file.wav", loop=True, mul=0.7)
pva = pyo.PVAnal(a, size=1024, overlaps=4, wintype=2)
pvs = pyo.PVSynth(pva).mix(2).out()
s.gui(locals())