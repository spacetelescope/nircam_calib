#########################
# Helpful Links
#########################

# Pipeline code
# https://github.com/spacetelescope/jwst/tree/master/jwst

# More on pytests
# https://docs.pytest.org/en/documentation-restructure/how-to/index.html




####################
Pytest dependencies for reporting and parallelization:
####################

conda install pytest-cov
conda install pytest-xdist
conda install pytest-html




#########################
# Notes
#########################
#
# ** unit test functions should not require or depend on input files
# ** if you think you need a reference file, you can use fake reference data
# ** all test function names must have a 'test_' prefix
# ** functions without the prefix won't be run as pytests (but other functions
#    can use them, e.g., setup_cube below)
# ** test names should be very descriptive for easy identification in the report




####################
Instructions for using tests with pipeline code:
####################

1.) Download git repository for the pipeline code.
2.) For each test you want to run, create a "tests" directory in the pipeline
    step directory (e.g., jwst/jwst/jump/tests) where you put your test code.
3.) From the pipeline step directory above the "tests" directory (e.g. "jump"
    directory above the "tests" directory), run pytest and it will check all
    code in the "jump" and "tests" directories.




####################
To run pytests:
####################

pytest tests/ --html=report.html --self-contained-html --cov=tests --cov-report=html


* The "html" keyword tells pytest to display the test results in an html file
  (it will be saved in the current directory).
* The "self-contained-html" keyword creates a self-contained html report instead
  of assets such as CSS and images getting stored separately by default to
  respect the Content Security Policy.
* The "cov" keyword tells pytest to show coverage for the routines in the listed
  directory.
* The "cov-report" keyword tells it to generate an html report that goes into a
  lower directory. The html for a given module is annotated to show what is and
  is not covered. To get the default formatting, you have to display the html
  report from a directory that has the needed javascript files (they end in .js)
  that are placed in the directory.




####################
To parallelize:
####################

pytest -n <cores>

with cores the number of machine cores to use for the tests (uses pytest-xdist)
