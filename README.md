# Musualizer

References:
- Music: https://www.youtube.com/@nu11_ft
- Raylib: https://www.raylib.com/
- FFT: https://www.youtube.com/watch?v=nmgFG7PUHfo
- Euler's Formula: https://en.wikipedia.org/wiki/Euler's_formula

## Audio Basics:

### Channels
- **Channels** = number of separate audio signals
- **1 channel** = Mono (single speaker)
- **2 channels** = Stereo (left + right speakers)

### Sample Size
- **Sample size** = bits used to represent one audio sample's amplitude
- Common values: 8-bit, 16-bit, 24-bit, **32-bit**
- Higher = better quality/precision

### Sample Rate
- How many samples per second (e.g., **44100 Hz** = CD quality)

### Frames
- **1 Frame = 1 sample per channel**
- For stereo (2 channels): **1 frame = 2 samples** (left + right)
- For mono (1 channel): 1 frame = 1 sample

Stereo Frame Layout (32-bit samples):
┌─────────────────────────────────────┐
│  Frame 0   │  Frame 1   │  Frame 2  │ ...
├─────────────────────────────────────┤
│ L0  │ R0   │ L1  │ R1   │ L2  │ R2  │ ...
│32bit│32bit │32bit│32bit │32bit│32bit│
└─────────────────────────────────────┘