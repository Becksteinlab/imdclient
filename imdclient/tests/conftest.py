"""
Global pytest fixtures
"""


# Command line arguments for 'test_manual.py'
def pytest_addoption(parser):
    parser.addoption(
        "--topol_arg",
        action="store",
    )
    parser.addoption(
        "--traj_arg",
        action="store",
    )
    parser.addoption("--first_frame_arg", action="store", type=int)


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    topol = metafunc.config.option.topol_arg
    traj = metafunc.config.option.traj_arg
    first_frame = metafunc.config.option.first_frame_arg

    if all(
        arg in metafunc.fixturenames
        for arg in ["topol_arg", "traj_arg", "first_frame_arg"]
    ):
        if topol is None or traj is None or first_frame is None:
            raise ValueError(
                "Must pass all three of '--topol_arg <path/to/topology>', "
                + "'--traj_arg <path/to/trajectory>', "
                + "'--first_frame_arg <first traj frame to compare to IMD>"
            )
        metafunc.parametrize("topol_arg", [topol])
        metafunc.parametrize("traj_arg", [traj])
        metafunc.parametrize("first_frame_arg", [first_frame])
