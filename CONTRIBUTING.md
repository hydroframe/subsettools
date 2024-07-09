# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

## Types of Contributions

### Report Bugs

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

### Running the test suite

We use [pytest](https://docs.pytest.org/) to run our tests. In general, a file
called `module.py` in our source code will have a corresponding `test_module.py`
in the test directory.

You can run the tests from the `subsettools` root directory with
```bash
pytest tests/
```

or just

```bash
pytest
```

To run the tests, you need to provide your HydroGEN PIN and email. To provide a
HydroGEN PIN for the tests, open the file subsettools/tests/conftest.py, 
uncomment lines 5 and 7 and provide your email and PIN.

After doing that, all the tests that check pure subsettools functionality should
pass. In addition, there are three test files (test_conus1_upper_verde.py,
test_conus1_upper_verde_spinup.py and test_multiple_hucs.py) that check
integration with ParFlow. These tests require an installation of ParFlow to run.
As ParFlow can be tricky to install, the easiest way would be to run them in our
ParFlow Docker container.

If you would like to run the ParFlow integration tests as well, you can:

1. Get our Docker container (https://hydroframesubsettools.readthedocs.io/en/latest/getting_started.html#using-docker)
2. Open a terminal inside the container and do:

```bash
$ pip install pytest
$ git clone https://github.com/hydroframe/subsettools.git
```

Then, you can edit the `conftest.py` file as described above and run `pytest`.

If you add a new feature, please make sure to include appropriate tests for that
feature in your pull request.

### Write Documentation

You can never have enough documentation! Please feel free to contribute to any
part of the documentation, such as the official docs, docstrings, or even
on the web in blog posts, articles, and such.

### Submit Feedback

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

## Get Started!

Ready to contribute? Here's how to set up `subsettools` for local development.

1. Download a copy of `subsettools` locally.
2. Install `subsettools` using `poetry`:

    ```console
    $ poetry install
    ```

3. Use `git` (or similar) to create a branch for local development and make your changes:

    ```console
    $ git checkout -b name-of-your-bugfix-or-feature
    ```

4. When you're done making changes, check that your changes conform to any code formatting requirements and pass any tests.

5. Commit your changes and open a pull request.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include additional tests if appropriate.
2. If the pull request adds functionality, the docs should be updated.
3. The pull request should work for all currently supported operating systems and versions of Python.

## Code of Conduct

Please note that the `subsettools` project is released with a
Code of Conduct. By contributing to this project you agree to abide by its terms.
