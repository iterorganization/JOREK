---
title: "Compress the Response Matrices"
nav_order: 7
parent: "Free Boundary Extension (STARWALL, CARIDDI)"
grand_parent: "Model Extensions"
layout: default
render_with_liquid: false
---

# Compress the Response Matrices

> **WARNING — THIS FEATURE HAS NOT BEEN MERGED.**
>
> Response-matrix compression is not available on the standard `develop` branch. The instructions on this page apply to the unmerged `feature/IMAS-4734-response-matrix-compression` branch and may differ from the final implementation. Do not expect these commands, files, or response-matrix formats to work in a normal JOREK checkout.

## Quick Overview

The implementation was developed on `feature/IMAS-4734-response-matrix-compression` and consists of two parts:

- The `compress_response` program, located in `jorek/vacuum/compression/`. Its compilation and use are described below.
- An adaptation of the `vacuum` module, located in `jorek/vacuum/`, which can automatically handle the compressed format.

An extended description of the method, objectives, and tests is available in the [response-matrix compression manuscript](https://arxiv.org/abs/2404.16546).

## Compiling and Running `compress_response`

The program currently has no dedicated input-file namelist. Its parameters must be edited directly in `jorek/vacuum/compression/mod_compression.f90`, as described in the next section.

`compress_response` uses ScaLAPACK. The `SCALAPACKLIB` and `SCALAPACKINCLUDE` variables must therefore be set in `Makefile.inc` before compilation. For example, the following configuration was used on MarconiFusion:

```makefile
##################
ifeq (1, ${COMPRESSION_TOOL})
  SCALAPACKLIB     = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
  SCALAPACKINCLUDE = -I${MKL_HOME}/include
endif
##################
```

Run the following command inside the `jorek` directory:

```bash
make -j8 compress_response
```

If compilation succeeds, it creates an executable named `compress_response`.

Run it from a submission script using:

```bash
mpirun -np ${NTASKS} ./compress_response
```

A `starwall-response.dat` file must be present in the run directory. With the SLURM scheduler, `NTASKS` can be replaced by `SLURM_NTASKS`. If the $n=0$ part of any matrix must remain uncompressed, a `boundary.txt` file must also be provided.

The program writes `starwall-response-compr.dat` with `file_version` set to:

- `7` if any matrix is factorized.
- `6` if all matrices use an aggregated format, including matrices re-aggregated after compression.

To use the result with JOREK, rename `starwall-response-compr.dat` to `starwall-response.dat`.

## Parameters in `mod_compression.f90`

Edit these parameters in `jorek/vacuum/compression/mod_compression.f90`. Default values are shown in bold.

| Parameter | Type | Possible values | Description |
|:---|:---|:---|:---|
| `compr_a_ye` | `logical` | `.true.` or **`.false.`** | Compress the `a_ye` response matrix. |
| `compr_a_ye_n0` | `logical` | `.true.` or **`.false.`** | Compress the $n=0$ part of `a_ye`. |
| `compr_a_ye_n0_sep` | `logical` | `.true.` or **`.false.`** | Compress each $n=0$ part of `a_ye` separately. |
| `compr_a_ey` | `logical` | `.true.` or **`.false.`** | Compress the `a_ey` response matrix. |
| `compr_a_ey_n0` | `logical` | `.true.` or **`.false.`** | Compress the $n=0$ part of `a_ey`. |
| `compr_a_ey_n0_sep` | `logical` | `.true.` or **`.false.`** | Compress each $n=0$ part of `a_ey` separately. |
| `compr_a_ee` | `logical` | `.true.` or **`.false.`** | Compress the `a_ee` response matrix. |
| `compr_s_ww` | `logical` | `.true.` or **`.false.`** | Compress the `s_ww` response matrix. |
| `compr_s_ww_inv` | `logical` | `.true.` or **`.false.`** | Compress the `s_ww_inv` response matrix. |
| `retained` | `real*8` | Used in `[0.d0, 1.d0]`; default **`-1.d0`** | Fraction of singular values retained after compression. |
| `retained_a_ye` | `real*8` | Used in `[0.d0, 1.d0]`; default **`-1.d0`** | Fraction of singular values retained for `a_ye`. |
| `retained_a_ey` | `real*8` | Used in `[0.d0, 1.d0]`; default **`-1.d0`** | Fraction of singular values retained for `a_ey`. |
| `retained_a_ee` | `real*8` | Used in `[0.d0, 1.d0]`; default **`-1.d0`** | Fraction of singular values retained for `a_ee`. |
| `retained_s_ww` | `real*8` | Used in `[0.d0, 1.d0]`; default **`-1.d0`** | Fraction of singular values retained for `s_ww`. |
| `retained_s_ww_inv` | `real*8` | Used in `[0.d0, 1.d0]`; default **`-1.d0`** | Fraction of singular values retained for `s_ww_inv`. |
| `debug` | `logical` | **`.true.`** or `.false.` | Print debug messages during compression. |
| `printsv` | `logical` | `.true.` or **`.false.`** | Write the singular values to a file. |
| `writeout` | `logical` | **`.true.`** or `.false.` | Write the output file. |
| `filevers` | `integer` | **`6`** or `7` | Select the output file version. |
| `CARIDDI` | `logical` | `.true.` or **`.false.`** | Read the input in `CARIDDI_mode`. |
| `reaggregate` | `logical` | **`.true.`** or `.false.` | Re-aggregate the factorized singular-value decomposition. |
| `verify_split` | `logical` | `.true.` or **`.false.`** | Check that the $n=0$ matrix splitting was performed correctly. |

## Vacuum-Module Notes and Limitations

### Notes

- The adapted vacuum module should automatically handle `starwall-response.dat` files with `file_version` equal to `6` or `7`.
- To run JOREK with compressed response matrices, rename `starwall-response-compr.dat` to `starwall-response.dat` and run JOREK normally.
- Several new variables are initialized in the `t_response_mat` data type defined in `jorek/vacuum/vacuum.f90`. See the comments in that file for details.

### Limitations

- The current implementation can handle compressed versions of `a_ye` and `a_ey` only.
- The formats of `a_ye` and `a_ey` must be identical: one matrix cannot be factorized while the other is aggregated.
