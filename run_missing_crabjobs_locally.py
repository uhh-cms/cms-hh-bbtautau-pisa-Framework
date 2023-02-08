import os
import sys
import json

from tqdm import tqdm
from argparse import ArgumentParser, RawDescriptionHelpFormatter

thisdir = os.path.realpath(os.path.dirname(__file__))

if not thisdir in sys.path:
    sys.path.append(thisdir)

import wlcg_dbs_interface as interface

verbosity = 0


def main(
    *args, 
    missing_files_json: str,
    **kwargs,
):
    # load dictionary with missing lfns
    missing_lfn_dict = dict()
    with open(missing_files_json) as f:
        missing_lfn_dict = json.load(f)
    
    pbar_samples = tqdm(missing_lfn_dict.keys())
    

def parse_arguments():
    description = """
    Script to run missing crab jobs locally and move the to the target
    (remote) site WLCG via gfal2.
    This script requires two things:
    - sample-config:        path to the sample config that maps a sample name to
                            the corresponding DAS key for the miniAOD
    - missing-files-json:   path to .json file containing the missing LFNs.

    The .json file with missing LFNs must have the following format:
    
    {
        "SAMPLE_NAME": {
            "missing": LIST_OF_LFNS,
            ...
        }
    }

    where 'SAMPLE_NAME' is also the key in the sample-config.
    This format is e.g. automatically produced with `check_crab_jobs.py'
    
    Please make sure to make gfal2 available to your python3 installation.
    On lxplus, you can do this with

    source "/cvmfs/grid.cern.ch/centos7-ui-160522/etc/profile.d/setup-c7-ui-python3-example.sh" ""

    Note that this only works on SLC7 machines at the moment, not on CentOS 8.
    """
    usage = "python %(prog)s [options] path/to/directories/containting/crab_base_dirs"


    parser = ArgumentParser(
        # usage=usage,
        description=description,
        formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument(
        "-w", "--wlcg-dir",
        help=" ".join("""
            path to your WLCG directory that is the final destination for your
            crab jobs. On T2_DESY, this would be your DCACHE (/pnfs) directory
        """.split()),
        metavar="PATH/TO/YOUR/WLCG/DIRECTORY",
        type=str,
        default=None,
        required=True,
        dest="wlcg_dir",
    )
    parser.add_argument(
        "--wlcg-prefix",
        help=" ".join(
            """
                Prefix to contact the WLCG directory at the remote site.
                Defaults to prefix for T2_DESY 
                (srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=)
            """.split()
        ),
        type=str,
        default="srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN=",
        dest="wlcg_prefix"
    )
    
    parser.add_argument(
        "--sample-config", "-c",
        help=" ".join(
            """
            path to sample config that contains the information about the
            original miniaod files (and thus the original name of the sample).
            Config must be the in yaml format!
            """.split()
        ),
        default=None,
        required=True,
        dest="sample_config",
        metavar="path/to/sample_config.yaml"
    )

    parser.add_argument(
        "-j", "--json",
        help=" ".join(
            """
                file that contains the missing lfns 
                (see structure in description)
            """.split()
        ),
        dest="missing_files_json",
        required=True,
        metavar="path/to/summary.json",
        type=str,
    )

    parser.add_argument(
        "-v", "--verbosity",
        type=int,
        default=0,
        help=" ".join(
            """
                control the verbosity of the output. Currently implemented levels
                0: just print number of files/outputs for each sample
                1: actually print the paths for the different files
            """.split()
        )
    )

    parser.add_argument(
        "-t", "--tmp-dir",
        help=" ".join(
            """
                directory to actually run the jobs in. Defaults to './tmp'
            """.split()
        ),
        type=str,
        dest="tmp_dir",
        metavar="path/to/temporary_directory",
        default="./tmp"
    )

    args = parser.parse_args()
    
    if not os.path.exists(args.sample_config):
        parser.error(f"file {args.sample_config} does not exist!")
    
    if not os.path.exists(args.missing_files_json):
        parser.error(f"file {args.missing_files_json} does not exist!")
    
    global verbosity
    verbosity = args.verbosity
    return args

if __name__ == '__main__':
    args = parse_arguments()
    main(**vars(args))