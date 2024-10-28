"""
Global pytest fixtures
"""


# Command line arguments for 'test_manual.py'
def pytest_addoption(parser):
    parser.addoption(
        "--topol_path_arg",
        action="store",
        default=None,
    )
    parser.addoption(
        "--traj_path_arg",
        action="store",
        default=None,
    )
    parser.addoption(
        "--first_frame_arg", action="store", type=int, default=None
    )


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    topol = metafunc.config.option.topol_path_arg
    traj = metafunc.config.option.traj_path_arg
    first_frame = metafunc.config.option.first_frame_arg

    if all(
        arg in metafunc.fixturenames
        for arg in ["topol_path_arg", "traj_path_arg", "first_frame_arg"]
    ):
        if topol is None or traj is None or first_frame is None:
            raise ValueError(
                "Must pass all three of '--topol_path_arg <path/to/topology>', "
                + "'--traj_path_arg <path/to/trajectory>', "
                + "'--first_frame_arg <first traj frame to compare to IMD>"
            )
        metafunc.parametrize("topol_path_arg", [topol])
        metafunc.parametrize("traj_path_arg", [traj])
        metafunc.parametrize("first_frame_arg", [first_frame])
