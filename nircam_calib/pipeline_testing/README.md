# Pipeline testing code

This repository contains code used to test the JWST calibration pipeline. Any contributed code should follow the `pytest` format. To read more about automating testing with `pytest`, visit the documentation pages: https://docs.pytest.org/en/documentation-restructure/how-to/index.html

There are several kinds of tests that will be structured differently described in the table below. 


| Test | Description |
| --- | --- |
| Unit | simple tests that don't require input files to check that pipeline functions do what they are expected to do |
| Regression | tests to see if pipeline steps produce expected outcome by comparing an input file with a truth file |
| Integration | end-to-end testing to check that all of the pipelines run as expected |
| Verification | functional testing to make sure pipeline runs for all types of data and modes without crashing (can be wrapped up into unit testing) |
| Validation P.1 | checks that pipeline meets requirements defined at start of mission planning (i.e., checks underlying algorithms) |
| Validation P.2 | determines accuracy and quality that pipeline achieves (to what extent selected algorithms meet error budget and how these vary for different types of data/science cases) |



If you would like to see the unit and regression tests that are already included in the pipeline, visit the pipeline GitHub repository: https://github.com/spacetelescope/jwst/tree/master/jwst. Unit tests are in the `tests` directories within the individual pipeline step subdirectories. Regression tests are in the `tests` directory listed in the `jwst` directory. 



## Using Python's `pytest` framework 

Below are some notes to help you get started with using `pytest`. To see an example script, take a look at `test_template_script.py` in this directory. 

#### Notes about tests

- unit test functions should not require or depend on input files (including reference files)
- all test function names must have a `test_` prefix (functions without the prefix won't be run as `pytests`, though other functions can use them)
- test names should be very descriptive for easy identification in the report


### Pytest dependencies for reporting and parallelization

To generate reports or take advantage of parallelization, you need to install the following plug-ins: 
- pytest-cov
- pytest-xdist
- pytest-html

Install them using conda or pip, e.g.: `conda install pytest-cov` or `pip install pytest-cov`. 

### To run `pytest`

This command depends on your specific directory structure, but this provides the main arguments: 

`$ pytest tests/ --html=report.html --self-contained-html --cov=tests --cov-report=html`

- `tests` directory contains the test scripts
- `html` keyword tells `pytest` to display the test results in an html file (it will be saved in the current directory).
- `self-contained-html` keyword creates a self-contained html report instead of storing assets such as CSS and images separately by default to respect the [Content Security Policy](https://developer.mozilla.org/en-US/docs/Web/HTTP/CSP).
- `cov` keyword tells `pytest` to show coverage for the routines in the listed directory.
- `cov-report` keyword tells `pytest` to generate an html report that goes into a lower directory. The html for a given module is annotated to show what is and is not covered. To get the default formatting, display the html report from a directory that has the javascript files (ending in `.js`) that are saved in the directory.


### To use the JWST pipeline `pytest` framework

1. Download the repository for the JWST pipeline code: https://github.com/spacetelescope/jwst
2. For each test you want to run, create a `tests` directory in the pipeline step directory (e.g., `jwst/jwst/jump/tests`) where you put your test code.
3. From the pipeline step directory above the `tests` directory (e.g., `jump` directory above the `tests` directory), run `pytest` command and it will check all the code in the pipeline step (`jump` and `tests`) directories.


### Parallelization
The `xdist` plug-in allows you to parallelize your tests. Do this by specifying the number of machine cores you want to use in place of `<cores>` below:

`$ pytest -n <cores>`
