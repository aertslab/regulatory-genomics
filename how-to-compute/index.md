# How to compute

Before starting the practicals you need a place to run the notebooks and a basic knowledge of the command line. The pages in this section explain how to connect to the teaching server, how to work in a Linux shell and how to install software in your own environment and use it from Jupyter.

The section consists of three parts, plus an extra page for students with a VSC account:

1. [Connecting to the teaching server](connect-to-teaching-server.md): logging in to the KU Leuven teaching server, either in the browser through JupyterHub or from a terminal with SSH (Linux, macOS, Windows).
2. [Shell tutorial](shelltutorial.ipynb): the basic shell commands for navigating folders and working with files (`mkdir`, `cd`, `ls`, `cp`, `mv`, `rm`, `gzip`, `head`, `wc`), redirecting and piping output, and searching and editing text files with `grep`, `tr`, `sed` and `awk`. The commands are meant to be typed in a terminal on the server, not run in the notebook. The tutorial is also available as a [PDF](shelltutorial.pdf).
3. [Environments and kernels](environment-and-kernel.ipynb): using the existing environments on the server, installing miniforge, creating your own environment with `mamba` and adding packages to it, and registering the environment as a Jupyter kernel so that a notebook can use it.
4. [Using the VSC](vsc.md) (only for students who already have a VSC account): the teaching server has CPUs only, so the VSC is the only way to use a GPU. During the sessions there are reserved CPU and GPU nodes and course credits, but you set up your own environments, kernels and data.

Every practical lists the software it needs and the environment and kernel it runs in on its topic page. If a notebook asks for a kernel that you do not have, part 3 shows how to make one.
