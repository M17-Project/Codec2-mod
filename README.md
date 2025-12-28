# Codec2-mod experimental fork
This repository contains a **minimal extraction of the Codec2's 3200 bps mode**, intended as a clean base for experimentation, optimization, and future research.

The initial goal of this work was to:
- Isolate only the functions and data structures actually used by the 3200 bps mode
- Achieve bit-exact encoder output compared to the reference `libcodec2`
- Remove unused code paths, state variables, and legacy scaffolding
- Establish a small codebase suitable for further development

Bit-exactness with the reference Codec2 encoder has been verified using identical input signals and byte-for-byte comparison of encoded frames.

## Motivation and goals

Once a minimal, bit-exact baseline is established, the planned next steps include:

- Removal of unused state and redundant computations
- Code optimizations
- Quantizer experiments (energy, pitch, LSPs)
- Exploration of improved excitation models
- Decoder-side enhancements that do not modify the Codec2 bitstream
- Suitability for embedded and low-power targets

This project is **not** intended to be a drop-in replacement for Codec2, but rather a controlled experimental platform derived from it.

## Branches

The `main` branch offers an encoder that is bitstream-compatible with vanilla Codec2 3200 bps mode.
The meaning, width, ordering, and allocation of all bit fields in the 3200 bps frame are all preserved, so bitstreams produced
by this encoder can be decoded by an unmodified Codec2 decoder (and vice versa).

The internal DSP implementation, floating-point operations, and decoded audio
signals are not required to be identical to the reference implementation.

Other branches may introduce experimental DSP changes (e.g. post-filters,
quantizers, or excitation models), potentially with a modified bitstream format.

## Important notice: derivative work

> [!NOTE]
> **This is a derivative work.**

This code is based heavily and directly on the Codec2 speech codec by venerable David Rowe, VK5DGR<sup>[1](https://github.com/drowe67) [2](https://www.qrz.com/db/VK5DGR)</sup> et al.

Original project:
- https://github.com/drowe67/codec2

Large portions of the code, algorithms, constants, and overall design originate from Codec2 and remain recognizably derived from it.  
All original credit for the Codec2 design, algorithms, and implementation belongs to David Rowe and the Codec2 contributors.

This repository exists to study, understand, and experimentally extend the Codec2 3200 bps mode.

## License

This project inherits the licensing requirements of Codec2.  
Please refer to the LICENSE file and original Codec2 license for details and ensure compliance when using or redistributing this code.

