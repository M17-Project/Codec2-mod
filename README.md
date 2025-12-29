# Codec2-mod experimental fork
This repository contains a minimal extraction of the Codec2's 3200 bps mode, intended as a clean base for experimentation, optimization, and future research.
Only the 3200 bps mode is supported - other Codec2 modes were intentionally excluded.

The initial goal of this work was to:
- Isolate only the functions and data structures actually used by the 3200 bps mode
- Achieve bit-exact encoder output compared to the reference `libcodec2`
- Remove unused code paths, state variables, and legacy scaffolding
- Establish a codebase that is as small and optimized as possible, suitable for further development

Bit-exactness with the reference Codec2 encoder has been verified using identical input signals and byte-for-byte comparison of encoded frames.
Bit-exactness refers to the encoded bitstream - decoded audio samples may differ from the reference implementation.
The removal of floating point maths was not a goal here.

## Motivation and goals

> [!NOTE]
> Bit-exact behavior compared to the vanilla Codec2 has been achieved.
> The `main` branch contains the refactored and optimized implementation while preserving bitstream compatibility.
> 
> The removal of floating-point math was not a goal of this work.

Planned next steps include:
- Quantizer experiments (energy, pitch, LSPs)
- Exploration of improved excitation models
- Decoder-side enhancements that do not modify the Codec2 bitstream (eg. the use of external neural networks)
- Fitness for embedded and low-power targets

This fork is intended to be a drop-in replacement for Codec2 when using the 3200 bps mode. It is also a controlled experimental platform derived from it.
The code structure and memory usage are intentionally designed to suit embedded systems.

## Current state and implemented changes

Compared to the reference Codec2 implementation, this fork already includes:
- Complete refactoring of the 3200 bps mode into a small codebase
- Removal of unused variables, modes, code paths, and legacy state not required for 3200 bps operation
- Elimination of all persistent dynamic memory allocation (no runtime `malloc`/`free`)
- Fully deterministic, fixed-size codec state suitable for static allocation
- Reduced overall memory footprint compared to the reference implementation
- Verified bitstream compatibility with the reference Codec2 encoder
- Initial ViSQOL evaluation shows a small MOS improvement in decoded audio compared to the reference implementation

These changes establish a stable and minimal baseline for further optimization and experimentation.

## Branches

The `main` branch offers an encoder that is bitstream-compatible with the reference Codec2 implementation.
The meaning, width, ordering, and allocation of all frame bit fields are all preserved. Bitstreams produced
by this encoder can be decoded by an unmodified Codec2 decoder (and vice versa).

The internal DSP implementation, floating-point operations, and decoded audio
signals are not required to be identical to the reference implementation.

Other branches may introduce experimental DSP changes (e.g. post-filters,
quantizers, or excitation models), potentially with a modified bitstream format.

## API differences vs. reference Codec2

This fork does not use the heap-allocated `codec2_create()` / `codec2_destroy()` API from the reference Codec2.
Instead, the codec state is explicitly owned by the caller and can be allocated statically.

### Reference Codec2 (libcodec2)

```c
struct CODEC2 *c2;

c2 = codec2_create(CODEC2_MODE_3200);

codec2_encode(c2, encoded, speech);
codec2_decode(c2, speech, encoded);

codec2_destroy(c2);
```

### Codec2-mod

```c
codec2_t c2;

codec2_init(&c2);

codec2_encode(&c2, encoded, speech);
codec2_decode(&c2, speech, encoded);
```

No destroy/free function is required, however `codec2_init()` has to be called before switching between encoder/decoder use.

## Important notice: derivative work

> [!NOTE]
> **This is a derivative work.**

This code is based heavily and directly on the Codec2 speech codec by venerable David Rowe, VK5DGR<sup>[1](https://github.com/drowe67) [2](https://www.qrz.com/db/VK5DGR)</sup> et al.

Original project:
- https://github.com/drowe67/codec2

Large portions of the code, algorithms, constants, and overall design originate from Codec2 and remain recognizably derived from it.  
All original credit for the Codec2 design, algorithms, and implementation belongs to David Rowe and the Codec2 contributors.

This repository exists to study, understand, optimize, and experimentally extend the Codec2 3200 bps mode.

## License

This project inherits the licensing requirements of Codec2.  
Please refer to the LICENSE file and original Codec2 license for details and ensure compliance when using or redistributing this code.

