# FAQ

## How do I install `subsettools` in a fresh virtual environment?

You can use conda to create an new virtual environment and pip to install subsettools:

```bash
conda create -n <YOUR_ENV_NAME> python=3.11 -y
conda activate <YOUR_ENV_NAME>
pip install subsettools
```

## What to do if I get a NumPy error when setting things up?

You may see the following error when setting up your environment:

  *A module that was compiled using NumPy 1.x cannot be run in*
  *NumPy 2.0.0 as it may crash. To support both 1.x and 2.x*
  *versions of NumPy, modules must be compiled with NumPy 2.0.*
  *Some module may need to rebuild instead e.g. with 'pybind11>=2.12'.*

  *If you are a user of the module, the easiest solution will be to*
  *downgrade to 'numpy<2' or try to upgrade the affected module.*
  *We expect that some modules will need time to support NumPy 2.*

NumPy is currently in transition from 1.x to 2.x. There are breaking changes between
the two versions, and they affect packages that are downstream from NumPy. There are
two solutions for this issue:

1. Downgrade the NumPy version to 1.26.4 with

```bash
pip install numpy==1.26.4
```
Make sure to restart your kernel after downgrading NumPy.

2. Make sure you have the latest version of `subsettools`, and upgrade if necessary:

You can check your installed version of `subsettools` with
```bash
pip show subsettools
```

You should compare the installed version with the latest published version on [PyPI](https://pypi.org/project/subsettools/).

If they differ, you can upgrade to the latest version of the package with
```bash
pip install subsettools --upgrade
```
Make sure to check that you now have the latest version installed.

Finally, it is possible that you have packages installed in your virtual environment with conflicting
requirements to `subsettools`. Some packages are still completing the migration to NumPy 2.x. If the two
solutions above fail, you can install `subsettools` in a fresh virtual environment as described in the
first FAQ above.
