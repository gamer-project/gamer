### What is in situ analysis?

In situ analysis is to analyze ongoing simulation data without outputting data to disk before conducting analysis.
It enables us to explore these ongoing simulation data through accessing memory directly using Python.

### Dependencies

- Install [libyt and yt_libyt](https://yt-project.github.io/libyt/HowToInstall.html#how-to-install): the in situ analysis library and its yt frontend.
- Install [yt](https://yt-project.org/): the core analysis tool.

### Compilation Options

libyt has the following compilation options:
- [[--libyt | Installation:-Option-List#--libyt]]
- [[--libyt_patchgroup | Installation:-Option-List#--libyt_patchgroup]]
- [[--libyt_interactive | Installation:-Option-List#--libyt_interactive]]
- [[--libyt_reload | Installation:-Option-List#--libyt_reload]]
- [[--libyt_jupyter | Installation:-Option-List#--libyt_jupyter]]

Must set `LIBYT_PATH` in the [[machine configuration file | Installation:-Machine-Configuration-File#1-library-paths]]
to the path containing the folders `include` and `lib` of `libyt`.

### Runtime Parameters

libyt has the following runtime parameters:
- `YT_SCRIPT`: The Python script name to be loaded. Do not include the file extension `.py`.
- `YT_VERBOSE`: Log level of libyt (0=off, 1=info, 2=warning, 3=debug).
- `YT_FIG_BASENAME`: Figure basename of the outputs from yt.

For example, `Input__Parameter` with the following lines will read the Python script `inline_script.py`,
set log level to information only, and set the figure basename to `MHD`.
```
YT_SCRIPT         inline_script
YT_VERBOSE        1
YT_FIG_BASENAME   MHD
```

### Quick Demos

1. [[Plummer | In-Situ-Python-Analysis:-Plummer]]

### FAQs

If you encounter the following OpenMPI error messages related to MPI remote memory access operation,
```
ompi_osc_ucx_win_attach: Assertion ......
```
try adding `OMPI_MCA_osc=sm,pt2pt` before `mpirun`:
```bash
OMPI_MCA_osc=sm,pt2pt mpirun ...
```

This is something libyt will fix in the near future.
