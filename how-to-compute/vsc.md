# Using the VSC (for students with a VSC account)

The teaching server only has CPUs. If you need a **GPU**, for example to train a model in the [s2f model training](../topics/s2f-model-training/index.md) practical, the only option is the Flemish Supercomputer Center (VSC). This page is for students who **already have a VSC account**. We cannot create a VSC account for just a single master course.

## Reserved nodes and credits

During the practical sessions, CPU and GPU nodes are reserved for this course. The compute time is paid from the course's credit account. Add these two options to every job you submit:

```bash
--reservation course_I0U29a_1      # the nodes reserved for the course (only during the sessions)
-A lp_edu_comparativegenomics      # the course's credit account
```

- **Batch job**: put the options in the header of your job script:

  ```bash
  #SBATCH --reservation=course_I0U29a_1
  #SBATCH -A lp_edu_comparativegenomics
  ```

- **Interactive job** (`srun`/`salloc`): add the options to the command, for example:

  ```bash
  srun --reservation=course_I0U29a_1 -A lp_edu_comparativegenomics --gpus-per-node=1 --pty bash -l
  ```

- **Jupyter through Open OnDemand** ([ondemand.hpc.kuleuven.be](https://ondemand.hpc.kuleuven.be)): when you start an interactive app (JupyterLab), fill in `lp_edu_comparativegenomics` as the account and `course_I0U29a_1` as the reservation. Ask for a GPU only if your notebook needs one.

The reservation only exists during the sessions. Outside the sessions, leave out `--reservation` and your job waits in the normal queue.

Use the credits fairly: stop interactive sessions and Jupyter servers when you are done, and only ask for a GPU when you are training a model.

## You are on your own here

On the VSC **none of the course environments, kernels or data are installed for you**. You have to set everything up yourself:

- **Environments and kernels**: install miniforge, create the environment listed on the topic page and register it as a Jupyter kernel, as explained in [Environments and kernels](environment-and-kernel.ipynb). Install miniforge and your environments in `$VSC_DATA`, not in your home folder: the home folder is too small.
- **Data**: the notebooks download their data or explain how to build the course data (see the topic pages). Keep large data in `$VSC_DATA` or `$VSC_SCRATCH`.
- **Paths**: the notebooks use paths relative to the notebook, so they run on the VSC as long as you keep the folder structure (`data/`, `results/` next to the notebook).

For questions about the VSC itself (logging in, quota, job scripts), see the [VSC documentation](https://docs.vscentrum.be).
