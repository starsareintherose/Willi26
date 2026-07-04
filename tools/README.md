# Tools

This directory contains small helper programs used for packaging of Willi26.

## Build release files

Cross compilation needs following dependencies installed:

* aarch64-linux-gnu-gcc: for cross compiling to aarch64 Linux

* osxcross: for cross compiling to macOS

* mingw-w64-gcc: for cross compiling to Windows

From the project root:

```sh
rustc tools/build_all.rs -O -o build_all
./build_all
```

It will generate binaries under `bin/`.


## Generate manual page

From the project root:

```sh
rustc tools/gen_man.rs -O -o gen_man
./gen_man
```

It will generate `man/willi26.1`
