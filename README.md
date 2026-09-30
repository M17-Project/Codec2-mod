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

Both implementations call the C library's `cosf()` during analysis, so the exact bitstream depends on its accuracy.
With an accurate `cosf()` (e.g. glibc), the encoded bitstreams are identical. With newlib-nano on a Cortex-M4,
the two differ in 1 of 150 frames of `hts1a`.

## Motivation and goals

> [!NOTE]
> Bit-exact behavior compared to the vanilla Codec2 has been achieved.
> The `main` branch contains the refactored and optimized implementation while preserving bitstream compatibility.
> 
> **The removal of floating-point math was not a goal of this work.**

Planned next steps include:
- Quantizer experiments (energy, pitch, LSPs)
- Exploration of improved excitation models
- Decoder-side enhancements that do not modify the Codec2 bitstream (eg. the use of external neural networks)
- Maintaining fitness for embedded and low-power targets

This fork is intended to be a drop-in replacement for Codec2 when using the 3200 bps mode. It is also a controlled experimental platform derived from it.
The code structure and memory usage are intentionally designed to suit embedded systems.

## Current state and implemented changes

Compared to the reference Codec2 implementation, this fork already includes:
- Complete refactoring of the 3200 bps mode into a small codebase
- Removal of unused variables, modes, code paths, and legacy state not required for 3200 bps operation
- Elimination of all persistent dynamic memory allocation (no runtime `malloc`/`free`)
- Fully deterministic, fixed-size codec state suitable for static allocation
- Much smaller flash footprint, no heap use, and a much smaller encoder stack than the reference implementation (see [Memory usage](#memory-usage))
- Verified bitstream compatibility with the reference Codec2 encoder

These changes establish a stable and minimal baseline for further optimization and experimentation.

## Speed comparison

STM32F411RE at 100 MHz, FPU enabled, *-Os* optimizations, newlib-nano, KISS FFT in all builds.
150 frames of `hts1a` (Codec2's standard test sample), timed per frame with the DWT cycle counter:

| Task                      | Reference Codec2 | Codec2-mod `main` | Codec2-mod `split-no-doubles` |
|---------------------------|------------------|-------------------|-------------------------------|
| Encoder init              | 18.07 ms¹        | 6.72 ms           | 3.08 ms                       |
| Decoder init              | 18.06 ms¹        | 5.00 ms           | 3.20 ms                       |
| Encode, avg / max         | 14.50 / 17.35 ms | 5.32 / 5.96 ms    | 4.17 / 4.55 ms                |
| Decode, avg / max         | 14.71 / 17.97 ms | 7.14 / 9.26 ms    | 5.47 / 6.61 ms                |
| Encode speedup (avg)      | 1×               | 2.7×              | 3.5×                          |
| Decode speedup (avg)      | 1×               | 2.1×              | 2.7×                          |
| CPU load, enc + dec (avg) | 146.0 %          | 62.3 %            | 48.1 %                        |

¹ The reference Codec2 has a single `codec2_create()` for both directions. It includes the heap allocation and freeing the previous instance with `codec2_destroy()`.

One frame lasts 20 ms. At 100 MHz, the reference Codec2 cannot encode and decode in real time simultaneously,
not even on average. Both Codec2-mod variants can, even in the worst-case frame (15.2 ms and 11.2 ms, respectively).

## Memory usage

Same target and build settings as above (`arm-none-eabi-gcc`, *-Os*, newlib-nano):

| Resource                          | Reference Codec2       | Codec2-mod `main` | Codec2-mod `split-no-doubles` |
|-----------------------------------|------------------------|-------------------|-------------------------------|
| Flash (code + constant data)²     | 167.2 KiB              | 23.8 KiB          | 19.5 KiB                      |
| Encoder state                     | 30.0 KiB³ (heap)       | 19.6 KiB (static) | 19.6 KiB (static)             |
| Decoder state                     | (same instance)³       | 16.2 KiB (static) | 16.2 KiB (static)             |
| Stack, encoder init               | 8.5 KiB                | 0.3 KiB           | 0.3 KiB                       |
| Stack, decoder init               | 8.5 KiB                | 0.1 KiB           | 0.1 KiB                       |
| Stack, `codec2_encode()`          | 14.1 KiB               | 1.0 KiB           | 1.0 KiB                       |
| Stack, `codec2_decode()`          | 14.4 KiB               | 5.1 KiB           | 5.1 KiB                       |
| Heap during encode/decode         | 0 B                    | 0 B               | 0 B                           |

**²** Taken from the linker map files. Includes the parts of libm, of the soft-float library and (for the reference Codec2) of the heap allocator that the codec pulls in. The reference Codec2 cannot link only the 3200 bps mode,
because `codec2_create()` references all modes.  
**³** One `codec2_create()` instance holds both the encoder and the decoder state. The figure includes the allocator's overhead.

For a full-duplex application, the total RAM (state + deepest stack) is about 44 KiB for the reference Codec2 and about 41 KiB for Codec2-mod. Codec2-mod needs no heap, and an encoder-only or decoder-only application
needs only the corresponding state.

## Branches

### `main`

The `main` branch offers an encoder that is bitstream-compatible with the reference Codec2 implementation.
The meaning, width, ordering, and allocation of all frame bit fields are all preserved. Bitstreams produced
by this encoder can be decoded by an unmodified Codec2 decoder (and vice versa).

The internal DSP implementation, floating-point operations, and decoded audio
signals are not required to be identical to the reference implementation.

### `split-no-doubles`
This branch gets rid of any double-precision arithmetic (usually by promotion).

> [!NOTE]
> Other branches may introduce experimental DSP changes (e.g. post-filters,
> quantizers, or excitation models), potentially with a modified bitstream format.

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

codec2_encode(&c2.encoder, encoded, speech);
codec2_decode(&c2.decoder, speech, encoded);
```

No destroy/free function is required. The encoder and decoder can be initialized and used independently.

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

