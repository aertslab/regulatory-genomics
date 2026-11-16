# Connecting to the Server

## A. Connect via browser directly in jupyter

### Go to https://jupyter.gbiomed.kuleuven.be/hub/

Fill in your r number and KUL (toledo) password

---

## B. Terminal

### Open or create your SSH config file (create a file by opening an unexisting file with nano)

```bash
nano ~/.ssh/config
```

### Paste the following into that file

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


### 1. linux


#### Download kmk (used for connecting to the server)

```bash
curl -s -o kmk https://repo.icts.kuleuven.be/artifactory/icts-p-lnx-generic-local/go/kmk/latest/kmk
```

#### running permissions

```bash
chmod 755 kmk
```

#### Run kmk

```bash
eval $(ssh-agent)
export KMK_USER="your r-number"
./kmk
```

#### sanity check

```bash
## this should not return an error or empty line
ssh-add -l
```

#### Connect

```bash
ssh gbw-s-bmw03
### alternative if above does not work:
ssh rXXXXXXX@gbw-s-bmw03
```

if prompted for password (KU Leuven password)

---

### 2. Terminal - MACOS

As in Linux, paste the config block above into `~/.ssh/config` (change it to your r-number!).

The easiest way forward is to download your kmk from the [KU Leuven SSH certificates page](https://admin.kuleuven.be/icts/services/ssh-cert/ssh-certificates-for-linux-mac#autotoc-item-autotoc-2).

```bash
## macOS with INTEL processor
curl -s -o kmk https://repo.icts.kuleuven.be/artifactory/icts-p-lnx-generic-local/go/kmk/latest/kmk-darwin-amd64
```

```bash
## macOS with M(ARM) processor
curl -s -o kmk https://repo.icts.kuleuven.be/artifactory/icts-p-lnx-generic-local/go/kmk/latest/kmk-darwin-arm64
```

#### Change permissions

```bash
chmod 755 kmk
```

#### Run kmk

```bash
eval $(ssh-agent)
export KMK_USER="your r-number"
./kmk
```

#### sanity check

```bash
## this should not return an error or empty line
ssh-add -l
```

#### Connect

```bash
ssh gbw-s-bmw03

### alternative if above does not work:
ssh rXXXXXXX@gbw-s-bmw03
```

if prompted for password (KU Leuven password)

---

### 3. Terminal - Windows

Send us an email if the above don't work.

---

Written by Gert Hulselmans (updated by Anton De Brabandere)