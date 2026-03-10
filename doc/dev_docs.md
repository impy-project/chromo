# `chromo` for developers

## Additional installation instructions

## From source

To install `chromo` from source, you will need to have a Python development environment set up, as well as C/C++ and Fortran compilers.

Begin by cloning the repository and checking out its submodules using the following command:

    git clone --recursive https://github.com/impy-project/chromo

Then, navigate to the cloned `chromo` directory:
    
    cd chromo

To install the package in editable mode (for developing the Python layer) and with verbose output, run the following command:

    pip install --no-build-isolation -v -e .[test,examples]

Please note that you will need to manually install the build environment in order for this command to succeed. You can check the `[build-system.requires]` key in [pyproject.toml](../pyproject.toml) to see which packages are required. Any warnings from the Fortran codes can be ignored.

Additionally, this command installs additional optional Python packages that are used in the tests and examples, but not required to run `chromo`.

## Known issues

- On OSX
    - You need to install `gcc` and `gfortran` with homebrew.
    - Apple introduced a bug in the Xcode Command Line Tools Version 14 which produces a linker error when compiling C++ code with `gcc`. Until this is fixed, the workaround is to downgrade to 13.4, use this link https://download.developer.apple.com/Developer_Tools/Command_Line_Tools_for_Xcode_13.4/Command_Line_Tools_for_Xcode_13.4.dmg and turn off automatic updates in the System Settings, because otherwise your Mac will upgrade to 14 again.

If you cannot fix the installation with these hints, please look into the subsection below which explains how to install in chromo in a verified docker environment. The docker environment has a properly set up environment verified by us, so that the installation is guaranteed to succeed.

## From source in Docker

This guide is intended for Linux and OSX users and assumes that you have a running Docker server. If you do not have Docker set up on your machine, please search online for instructions on how to do so.

To begin, clone the `chromo` repository and navigate to the cloned directory:

    git clone --recursive https://github.com/impy-project/chromo
    cd chromo

Next, download the Linux image for x86_64 by running the following command:

    docker pull quay.io/pypa/manylinux2014_x86_64

If you are running aarch64 or a VM on Apple Silicon, use the following image instead and adjust the end of the next command accordingly:

    docker pull quay.io/pypa/manylinux2014_aarch64

After downloading the image, create a Docker instance and bind the `chromo` directory to it by running the following command:

    docker run --rm -d -it --name chromo -v "$(pwd)":/app quay.io/pypa/manylinux2014_x86_64

This command enters your Docker instance and changes the working directory to /app. From there, select your desired Python version (e.g., 3.11) and enter a virtual environment:

    cd /app
    python3.11 -m venv venv
    source venv/bin/activate

Finally, install chromo and its dependencies using the following command:

    pip install --prefer-binary -v -e .

Now, you can use `chromo` within the Docker instance. If you are using Linux, you can also create a wheel inside Docker and install it on your host machine by following these steps:

    # inside Docker
    pip install wheel
    python setup.py bdist_wheel

    # exit Docker with ctrl+D
    pip install dist/*.whl

Please note that this will only work if you are using the same Python version inside and outside of Docker.

## Running tests

Some notes regarding tests.

- Tests are run in parallel by default with `pytest-xdist`. To disable this, use the option `-n 0`.
- The test `test_generators` takes a long time. You can skip it with the option `-k "not test_generators"`.
- Tests which run a model do so in a separate process, because most models can only be instantiated once. This prevents using `--pdb` to start the debugger at the point of failure. You can prefix the `pytest` call like this `DEBUG=10 python -m pytest ...` to run the model in the current process. This will only work once for each model and lead to failures afterwards.