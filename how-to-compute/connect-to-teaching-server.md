# Connecting to the Server

## Linux

### 1. Open or create your SSH config file (create a file by opening an unexisting file with nano)

```bash
nano ~/.ssh/config
```

### 2. Paste the following into that file

```ssh-config
Host gbw-s-* !jump-s.icts.kuleuven.be
    User rXXXXXXX
    HostName %h.luna.kuleuven.be
    ProxyJump rXXXXXXX@jump-s.icts.kuleuven.be

Host gbw-s-bmw0*
    PubkeyAuthentication no
    PreferredAuthentications password
```

> **Note:** Replace `rXXXXXXX` with your own r-number.
> **Note:** VSC users should already have something in their ssh config

### 3. Download kmk (used for connecting to the server)

```bash
wget -O kmk https://repo.icts.kuleuven.be/artifactory/icts-p-lnx-generic-local/go/kmk/latest/kmk
```

### 4. Set up

```bash
eval $(ssh-agent)
export KMK_USER="rXXXXXXX"
```

### 5. Give kmk running permissions

```bash
chmod u+x kmk
```

### 6. Run kmk

```bash
./kmk
```

### 7. Connect

```bash
ssh gbw-s-bmw03
```

---

## macOS

As in Linux, paste the config block above into `~/.ssh/config` (change it to your r-number!).

The easiest way forward is to download your kmk from the [KU Leuven SSH certificates page](https://admin.kuleuven.be/icts/services/ssh-cert/ssh-certificates-for-linux-mac#autotoc-item-autotoc-2).

```bash
wget https://admin.kuleuven.be/icts/services/ssh-cert/kmk-darwin-arm64 # MACOS M1/2/3/4
wget https://admin.kuleuven.be/icts/services/ssh-cert/kmk-darwin-amd64 # MACOS Intel
wget https://admin.kuleuven.be/icts/services/ssh-cert/kmk              # Linux
```

### 1. Change permissions

```bash
chmod +x kmk-darwin-arm64
```

### 2. Run kmk

```bash
./kmk-darwin-arm64
```

### 3. Connect

```bash
ssh gbw-s-bmw03
```

---

## Windows

Send us an email if the above don't work.

---

Written by Gert Hulselmans