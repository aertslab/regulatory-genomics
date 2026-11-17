# Connecting to the Server

## A. Connect via browser directly in jupyter

### Go to https://jupyter.gbiomed.kuleuven.be/hub/

Fill in your r number and KUL (toledo) password

---

## B. Terminal

The server can only be reached through a KU Leuven jump server, with a certificate that you get after logging in with 2FA. The tool `kmk` gets this certificate for you. You only have to do steps 1 and 2 once.

> **Note:** Replace `rXXXXXXX` with your own r-number everywhere below.

### 1. Download kmk (once)

Run **one** of these lines, depending on your computer:

```bash
## Linux
curl -s -o kmk https://repo.icts.kuleuven.be/artifactory/icts-p-lnx-generic-local/go/kmk/latest/kmk

## macOS with INTEL processor
curl -s -o kmk https://repo.icts.kuleuven.be/artifactory/icts-p-lnx-generic-local/go/kmk/latest/kmk-darwin-amd64

## macOS with M(ARM) processor
curl -s -o kmk https://repo.icts.kuleuven.be/artifactory/icts-p-lnx-generic-local/go/kmk/latest/kmk-darwin-arm64
```

Give it running permissions:

```bash
chmod 755 kmk
```

### 2. Add the server to your SSH config (once)

Open or create your SSH config file (create a file by opening an unexisting file with nano):

```bash
nano ~/.ssh/config
```

Paste the following into that file:

```ssh-config
Host gbw-s-* !jump-s.icts.kuleuven.be
    User rXXXXXXX
    HostName %h.luna.kuleuven.be
    ProxyJump rXXXXXXX@jump-s.icts.kuleuven.be

Host gbw-s-bmw0*
    PubkeyAuthentication no
    PreferredAuthentications password
```

> **Note:** VSC users should already have something in their ssh config

### 3. Get a certificate (once a day)

```bash
eval $(ssh-agent)
export KMK_USER="rXXXXXXX"
./kmk
```

A web page opens: log in with your KU Leuven account and 2FA.

Sanity check, this should show a line ending in `(RSA-CERT)`, not an error or an empty line:

```bash
ssh-add -l
```

The certificate is valid for one day. The next day, or when you open a new terminal, run these lines again.

### 4. Connect

```bash
ssh gbw-s-bmw03
### alternative if above does not work:
ssh rXXXXXXX@gbw-s-bmw03
```

If prompted for a password, use your KU Leuven password.

If it does not connect from home, connect to the KU Leuven VPN first and try again.

---

### Windows

Install the certificate tool `certagent` from the Software Center (or download [certagent.exe](https://repo.icts.kuleuven.be/artifactory/icts-p-lnx-generic-local/go/certagent/latest/certagent.exe)), run it and right-click it to generate your certificate. Then connect with a program that supports jump servers (Bitvise, or PuTTY), using `jump-s.icts.kuleuven.be` as jump server and `gbw-s-bmw03.luna.kuleuven.be` as server.

Send us an email if this does not work.

---

Written by Gert Hulselmans (updated by Anton De Brabandere)
