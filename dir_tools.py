import os
import sys
import re


def get_all_paths():
    run_dir = ""

    try:
        run_dir = sys.argv[1]
        if not os.path.isdir(run_dir):
            raise OSError("Not a valid path.")
    except IndexError:
        print("""Please give the path to the run directory as a command line
               argument.""")
        sys.exit(0)
    except OSError as err:
        print(err.args)
        sys.exit(0)

    # adding forward slash at end
    if run_dir[-1] is not "/":
        run_dir = run_dir + "/"

    # get run name i.e. last dir
    run_name = ""
    run_re = re.match(r"[\S*\s*]*/([\S*\s*]*)/\Z", run_dir)

    run_name = run_re.group(1)
    root = os.getcwd() + "/"
    return ({"root": root,
             "run_name": run_name,
             "run_dir": run_dir})
