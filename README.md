# Musializer

A real-time audio frequency spectrum visualizer using FFT (Fast Fourier Transform).

## Features

- Real-time FFT-based frequency spectrum visualization
- 256 frequency bins
- Play/pause with spacebar
- Dark theme UI

## Building

```bash
./build.sh
```

## Usage

```bash
./build/musializer <audio-file>
```

Example:
```bash
./build/musializer music/song.mp3
```

## Controls

| Key | Action |
|-----|--------|
| Space | Play/Pause |
| ESC | Exit |

## Project Structure

```
musializer/
├── src/
│   ├── musializer.c    # Main application
│   └── fft.h           # FFT/DFT implementation (header-only)
├── build.sh            # Build script
└── README.md
```

## Dependencies

- [raylib](https://www.raylib.com/) - Graphics and audio library
- GLFW
- C99 compiler (clang/gcc)

## References

- Music: https://www.youtube.com/@nu11_ft
- Raylib: https://www.raylib.com/
- FFT: https://www.youtube.com/watch?v=nmgFG7PUHfo
- Euler's Formula: https://en.wikipedia.org/wiki/Euler's_formula

---

## Audio Basics

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

```
Stereo Frame Layout (32-bit samples):
┌─────────────────────────────────────┐
│  Frame 0   │  Frame 1   │  Frame 2  │ ...
├─────────────────────────────────────┤
│ L0  │ R0   │ L1  │ R1   │ L2  │ R2  │ ...
│32bit│32bit │32bit│32bit │32bit│32bit│
└─────────────────────────────────────┘
```

---

## FFT Theory

### What is FFT?

The **Fast Fourier Transform** decomposes a signal into its constituent frequencies.

- **Input**: Time-domain signal (amplitude vs time)
- **Output**: Frequency-domain spectrum (amplitude vs frequency)

### DFT vs FFT

| Aspect | DFT | FFT |
|--------|-----|-----|
| Complexity | O(N²) | O(N log N) |
| N=1024 | 1,048,576 ops | 10,240 ops |
| N=4096 | 16,777,216 ops | 49,152 ops |

### Key Concepts

- **Frequency Resolution** = Sample Rate / N
  - 44100 Hz / 256 = ~172 Hz per bin
- **Nyquist Frequency** = Sample Rate / 2
  - Maximum detectable frequency (22050 Hz for CD audio)
