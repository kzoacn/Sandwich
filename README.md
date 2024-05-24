# Sandwich

This code is modified from [FAEST](https://github.com/faest-sign/faest-ref)

Please install openssl, otherwise the execution will be very slow.

If you are using Ubuntu OS, run the following commands to install some packages.

```
sudo apt install meson libssl-dev libboost-dev
```

## Run

```
meson setup build
cd build
meson configure -Dbenchmarks=enabled
cd ..
bash bench.sh
```
