import numpy as np
import pytest

# @pytest.mark.hookwrapper
# def pytest_runtest_makereport(item, call):
#     pytest_html = item.config.pluginmanager.getplugin('html')
#     outcome = yield
#     report = outcome.get_result()
#     extra = getattr(report, 'extra', [])
#     if report.when == 'call':
#         # always add url to report
#         extra.append(pytest_html.extras.url('http://www.example.com/'))
#         xfail = hasattr(report, 'wasxfail')
#         if (report.skipped and xfail) or (report.failed and not xfail):
#             # only add additional html on failure
#             extra.append(pytest_html.extras.html('<div>Additional HTML</div>'))
#         report.extra = extra

def pytest_generate_tests(metafunc):
    if 'cases' in metafunc.funcargnames:
        metafunc.parametrize("cases", [
            ("NRCNRCA1-DARK-60012216201_1_481_SE_2016-01-02T02h34m28_uncal_sliced_dq_init_saturation_superbias.fits")
        ])

    if 'smoothing_lengths' in metafunc.funcargnames:
        metafunc.parametrize("smoothing_lengths", [
            (11)
        ])

    # absolute tolerance, relative tolerance
    # for testing abs(manual - pipeline) <= abs_tol + rel_tol * abs(pipeline)
    # if equation == True, test is passed
    if 'tolerances' in metafunc.funcargnames:
        metafunc.parametrize("tolerances", [
            ([1e-2,1e-2])
        ])

    if 'sigmas' in metafunc.funcargnames:
        metafunc.parametrize("sigmas", [
            (3)
        ])

    if 'side_gains' in metafunc.funcargnames:
        metafunc.parametrize("side_gains", [
            (1.)
        ])

    # # absolute tolerance, relative tolerance
    # # for testing abs(manual - pipeline) <= abs_tol + rel_tol * abs(pipeline)
    # # if equation == True, test is passed
    # if 'tolerances' in metafunc.funcargnames:
    #     metafunc.parametrize("tolerances", [
    #         ([1e-2,1e-3]),
    #         ([1e-2,1e-2])
    #     ])

    # if 'sigmas' in metafunc.funcargnames:
    #     metafunc.parametrize("sigmas", [
    #         (2),
    #         (3),
    #         (4)
    #     ])

    # if 'side_gains' in metafunc.funcargnames:
    #     metafunc.parametrize("side_gains", [
    #         (1.),
    #         (1.5),
    #         (2.)
    #     ])
