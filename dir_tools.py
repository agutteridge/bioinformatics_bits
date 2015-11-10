import os
import sys
import re


def get_all_paths(args_list):
    run_dir = ""
    run_name = ""

    # Requires path to run directory to be supplied as a command line argument.
    try:
        run_dir = args_list[1]

        if not os.path.isdir(run_dir):
            raise OSError("Not a valid path.")

        # adding forward slash at end
        if run_dir[-1] is not "/":
            run_dir = run_dir + "/"

        # get run name i.e. last dir
        # Assumes Illumina numeric ID is at least 8 digits long
        run_re = re.match(r"([\S*\s*]*)-\d\d\d\d\d\d\d\d/+\Z", run_dir)

        if run_re:
            run_name = run_re.group(1)
        else:
            raise AttributeError("Not a run directory.")
    except IndexError:
        print("Please give the path to the run directory as a command line" +
              " argument.")
        sys.exit(0)
    except OSError as oe:
        for e in oe.args:
            print(e)
        sys.exit(0)
    except AttributeError as ae:
        for e in ae.args:
            print(e)
        sys.exit(0)

    root = os.getcwd() + "/"
    return ({"root": root,
             "run_name": run_name,
             "run_dir": run_dir})
