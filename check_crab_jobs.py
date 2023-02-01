import os
import sys
import json

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from tqdm import tqdm
from itertools import chain
from collections.abc import Iterable

def check_job_outputs(
    wlcg_dir: str,
    wlcg_prefix: str,
    timestamp: str,
    collector_set: set[str],
    input_map: dict[str, list[str]],
    job_details: dict[str, dict],
    state: str="failed"
) -> None:
    """Function to collect information about jobs in *job_details*.
    First, all job ids with state *state* are retrieved from *job_details*.
    Then, the *collector_set* is filled with information depending on the 
    *state*.
    If *state* is 'failed', the *collector_set* is filled with paths to files
    that correspond to failed jobs and thus should not exist.
    If *state* is 'finished', the *collector_set* is filled with lfns that
    were successfully processed. In this case, an additional check whether
    a lfn is already marked as done is performed, and raises an error if
    a marked lfn is supposed to be added again.

    TODO: implement file readout with gfal

    Args:
        wlcg_dir (str): path for WLCG directory
        wlcg_prefix (str): prefix to contact the WLCG directory via gfal
        timestamp (str):    time stamp of the job to find the correct path
                            on the (remote) WLCG site
        collector_set (set):    set to be filled with information, depending on 
                                *state* (see description above)
        input_map (dict): Dictionary of format {job_id: list_of_lfns}
        job_details (dict): Dictionary containing the status of the jobs of
                            format {job_id: ADDITIONAL_INFORMATION}.
                            *ADDITIONAL_INFORMATION* is a dict that must contain
                            the keyword 'State'.
        state (str, optional):  State to select the relevant job ids.
                                Must be either "failed" or "finished".
                                Defaults to "failed".

    Raises:
        ValueError: If a lfn is already marked as done but is associated with
                    a done job again, the ValueError is raised.
    """
    relevant_ids = filter(
        lambda x: job_details[x]["State"] == state,
        job_details
    )
    #TODO implement readout via gfal


    # for state "failed", collect output files that should not be there
    if state == "failed":
        pass
    # if state is finished, safe the done lfns (if the output of the job is also 
    # available)
    elif state == "finished":
        # first check if a lfn is already marked as done - this should not happen
        lfns = set(chain.from_iterable([input_map[x] for x in relevant_ids]))
        overlap = collector_set.intersection(lfns)
        if len(overlap) != 0:
            overlap_string = "\n".join(overlap)
            raise ValueError(f"""
            The following lfns were already marked as done:
            {overlap_string}

            This should not happen!
            """)
        
        collector_set.update(lfns)

def get_job_inputs(crab_dir: str, job_input_file: str="job_input_files.json"):
    """Load job input file. This .json file contains the mapping of the form
    {
        'job_id': [
            'lfn_1',
            'lfn_2',
            ...
        ],
        ...
    }

    Args:
        crab_dir (str): current crab directory
        job_input_file (str, optional): file name in *crab_dir* containing the
                                        mapping. 
                                        Defaults to "job_input_files.json".
    
    Raises:
        ValueError: If the file containing the file mapping does not exist

    Returns:
        dict: file mapping of format described above
    """    

    # first, build file path
    path = os.path.join(crab_dir, job_input_file)

    # if the file does not exist, raise an error
    if not os.path.exists(path):
        raise ValueError(f"Could not load input mapping from file '{path}'")

    # load json file
    with open(path) as f:
        input_map = json.load(f)

    return input_map

def get_status(sample_dir: str, status_file: str, crab_dir: str) -> dict[str, dict]:
    """Function to load the status from a .json file in *sample_dir*.
    If the file 'sample_dir/status_file.json' does not exist, the script falls
    back to loading 'status.json'

    Args:
        sample_dir (str): path to directory containing the status json file
        status_file (str): name of the .json file containing the job stati
        crab_dir (str): path to current crab base directory - currently not used

    Raises:
        NotImplementedError:    If the job stati cannot be loaded from a json
                                file they must be obtained with `crab status`
                                directoy. This is not implemented at the moment


    Returns:
        dict: Dictionary with status information given by crab
    """    

    # build path to the file containing the status information
    status_path = os.path.join(sample_dir, f"{status_file}.json")

    # also build a fallback in case the current file does not exist
    backup_status_path = os.path.join(sample_dir, "status.json")

    # now check if the file in either *status_path* or the backup exists
    if os.path.exists(status_path):
        with open(status_path) as f:
            status = json.load(f)
    elif os.path.exists(backup_status_path):
        with open(backup_status_path) as f:
            status = json.load(f)
    else:
        # if they do not exist, raise an error
        raise NotImplementedError("Obtaining the status from crab not implemented yet!")
    
    return status

def check_status(status: dict[str, dict], crab_dir: str) -> None:
    """Function to ensure that the *status* dictionary indeed belongs to the
    current crab base directory *crab_dir*. To this end, the 'project_dir'
    in the status is compared to the absolute path of *crab_dir*.

    Args:
        status (dict): Dictionary containing status information about the jobs
        crab_dir (str): path to current crab base directory

    Raises:
        ValueError: if the 'project_dir' does not match the absolute path of
                    the current crab base directory, an error is raised
    """    
    # load project dir from status
    status_project_dir = status.get("project_dir", None)

    # build absolute path for *crab_dir*
    abs_crab_dir = os.path.abspath(crab_dir)

    # do the comparison
    if not status_project_dir or not status_project_dir == abs_crab_dir:
        raise ValueError("Project dir in status file does not match current crab dir under scutiny!")

def check_crab_directory(
    sample_dir: str,
    sample_name: str,
    suffix: str,
    status_file: str,
    known_lfns: set[str],
    done_lfns: set[str],
    failed_job_outputs: set[str],
    pbar: Iterable,
    wlcg_dir: str,
    wlcg_prefix: str,
    job_input_file: str="job_input_files.json"
) -> None:
    """Function to check a specific crab base directory in *sample_dir*.
    First, the name of the crab base directory (*crab_dir*) is built from
    *sample_dir*, *sample_name* and *suffix*.
    Then, the input file list and the status of the jobs of *crab_dir* is loaded.
    For book-keeping purposes, a set of known lfns *known_lfns* is also required.
    If this set is empty (which should happen for the first crab base directory
    that (should) consider all lfns), fill this set.
    If there are unknown lfns in the other (recovery) jobs, raise an error.
    Additionally, obtain the time stamp of the current crab job to correctly
    load/check the output of the jobs on the WLCG site.
    Finally, the relevant information about failed and finished jobs is obtained.

    Args:
        sample_dir (str): path to the directory containing the crab base directories
        sample_name (str): Name of the current sample (as it is in sample_dir).
                            Used to build name of crab directory via 'crab_*sample_name*'
        suffix (str):   suffix to use for the current crab directory, such that 
                        the final name of the *crab_dir* is 
                        'crab_*sample_name*_*suffix*'
        status_file (str):  name of the status file containing the information 
                            for the current crab base directory *crab_dir* 
        known_lfns (set[str]): Set of already known lfns
        done_lfns (set[str]): Set of lfns corresponding to successful jobs
        failed_job_outputs (set[str]):  Set of output files originating from 
                                        failed jobs. This should not happen,
                                        thus we track them
        pbar (Iterable): Status bar for the iteration in the main function.
                                Used to correctly display the current *crab_dir*.
        wlcg_dir (str): Name of the WLCG Directory containing the outputs 
                        of the jobs.
        wlcg_prefix (str): Prefix to contact the WLCG Directory with gfal
        job_input_file (str, optional): Name of the file containing the
                                        mapping of job_id -> input file(s) for
                                        a given *crab_dir*. 
                                        Defaults to "job_input_files.json".

    Raises:
        ValueError: If previously unkown lfns are encountered
        ValueError: if the time stamp of a given crab job cannot be obtained
    """    
    # build name of current crab base directory
    if not suffix == "" and not suffix.startswith("_"):
        suffix = "_"+suffix
    crab_dirname = f"crab_{sample_name}"+suffix
    crab_dir = os.path.join(sample_dir, crab_dirname)
    if not os.path.exists(crab_dir):
        print(f"Directory {crab_dir} does not exist, will stop looking here")
        return
    pbar.set_description(f"Checking directory {crab_dir}")
    
    # load the input file mapping of the form 'job_ids' -> list of lfns
    input_map = get_job_inputs(crab_dir=crab_dir, job_input_file=job_input_file)
    
    # perform sanity checks
    # first load a flat list of lfns
    flat_lfns = set(chain.from_iterable(input_map.values()))

    # if we don't know any lfns yet, use this set as a baseline
    if len(known_lfns) == 0:
        known_lfns.update(flat_lfns)
    
    # check if all lfns are known at this point
    unknown_lfns = flat_lfns.difference(known_lfns)
    if len(unknown_lfns) != 0:
        # if we enter here, lfns appeared that were previously unknown
        # this shouldn't be possible (unless maybe due to TAPE_RECALLS)
        # currently being checked
        unknown_lfns_string = "\n".join(unknown_lfns)
        from IPython import embed; embed()
        raise ValueError(f"""
            Following lfns are not known in '{crab_dir}':
            {unknown_lfns_string}

            This should not happen!
            """)

    # load the dictionary containing the job stati
    status = get_status(
        sample_dir=sample_dir,
        status_file=status_file,
        crab_dir=crab_dir
    )
    # sanity check whether we indeed loaded the correct status file
    check_status(status=status, crab_dir=crab_dir)

    # get general information about jobs
    n_jobs = status.get("n_jobs_total", 0)

    # load the information about the specific jobs that form the crab job
    job_details = status.get("details", dict())

    # load time stamp
    time_stamp = status.get("task_name", None)
    if not time_stamp:
        raise ValueError("Could not retrieve time stamp from status json!")
    time_stamp = time_stamp.split(":")[0]

    # load information about failed jobs
    check_job_outputs(
        wlcg_dir=wlcg_dir,
        wlcg_prefix=wlcg_prefix,
        timestamp=time_stamp,
        collector_set=failed_job_outputs,
        input_map=input_map,
        job_details=job_details,
        state="failed",
    )

    # load information about failed jobs
    check_job_outputs(
        wlcg_dir=wlcg_dir,
        wlcg_prefix=wlcg_prefix,
        timestamp=time_stamp,
        collector_set=done_lfns,
        input_map=input_map,
        job_details=job_details,
        state="finished",
    )


def main(*args, **kwargs):
    """main function. Load information provided by the ArgumentParser. Loops
    Thorugh the sample directories provided as *sample_dirs* and the *suffices*
    to check the individual crab base directories.
    Finally, check if any lfns are unaccounted for in the list of finished jobs.
    """  
    # load the information from the argument parser  
    sample_dirs = kwargs.get("sample_dirs", [])
    wlcgs_dir = kwargs.get("wlcg_dir", None)
    wlcg_prefix = kwargs.get("wlcg_prefix", "")
    suffices = kwargs.get("suffix", [])
    status_files = kwargs.get("status_files", [])

    missing_lfns = dict()
    # loop through the sample directories containing the crab base directories
    pbar_sampledirs = tqdm(sample_dirs)
    for sample_dir in pbar_sampledirs:
        # if a sample dir does not exist, no need to check it
        if not os.path.exists(sample_dir):
            print(f"Directory {sample_dir} does not exist, skipping!")
            continue
        sample_dir = sample_dir.strip(os.path.sep)

        # from IPython import embed; embed()
        # extract the sample name from the sample directory
        sample_name=os.path.basename(sample_dir)
        pbar_sampledirs.set_description(f"Checking sample {sample_name}")

        # set up the sets to keep track of the lfns
        known_lfns=set()    # full set of lfns for this sample
        done_lfns=set()     # set of lfns processed by successful jobs
        
        # set of **outputs** from failed jobs, which shouldn't happen
        failed_job_outputs=set()    

        # loop through suffices to load the respective crab base directories
        pbar_suffix = tqdm(zip(suffices, status_files))
        for suffix, status_file in pbar_suffix:
            check_crab_directory(
                sample_dir=sample_dir,
                sample_name=sample_name,
                suffix=suffix,
                status_file=status_file,
                known_lfns=known_lfns,
                done_lfns=done_lfns,
                failed_job_outputs=failed_job_outputs,
                pbar=pbar_suffix,
                wlcg_dir=wlcgs_dir,
                wlcg_prefix=wlcg_prefix,
            )

        # in the end, all LFNs should be accounted for
        unprocessed_lfns = known_lfns.symmetric_difference(done_lfns)
        if len(unprocessed_lfns) != 0:
            missing_lfns[sample_dir] = unprocessed_lfns
    
    # do some final sanity check: if we found missing lfns, print them here
    # so the user can do something
    if len(missing_lfns) != 0:
        print("found missing lfns for the following samples")
        print(missing_lfns)
    else:
        print("everything ok!")


def parse_arguments():
    description = """
    Small script to cross check crab jobs. This script checks the following:
    - Was every LFN entry in DAS processed?
    - Are there outputs for all jobs with status "DONE"?
    - Are there outputs for jobs with status "FAILED"?
    - Is there an output for every LFN?

    The script expects the following directory structure:
    - SAMPLE_NAME
    |   ---crab_SAMPLE_NAME
    |   ---crab_SAMPLE_NAME_recovery_1
    |   ---crab_SAMPLE_NAME_recovery_2
    |   ---...
    - NEXT_SAMPLE
    ....

    This structure is created by the crabOverseer by default.
    You can specify which suffixes (e.g. *recovery_1*) are to be checked, 
    see options below. Note that the order you specify in this option is also
    the order in which the directories are checked.
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
        default="srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN="
    )
    parser.add_argument(
        "-s", "--suffix",
        help=" ".join("""
            specify the suffices you would like for the check of the crab
            source directory. Can be a list of suffices, e.g. 
            `-s "" recovery_1 recovery_2`. Defaults to 
            ["", "recovery_1", "recovery_2`"]
        """.split()),
        default=None,
        nargs="+",
    )
    parser.add_argument(
        "--status-files",
        help=" ".join(
            """
            List of status files containing job information.
            Entries in this list are zipped to the list of suffixes (option `-s`)
            Therefore, the length of this list must be the same as list
            obtained from option `--suffix`! 
            Defaults to ["status_0", "status_1", "status"]
            """.split()
        ),
        default=None,
        nargs="+",
        dest="status_files"
    )

    parser.add_argument(
        "sample_dirs",
        help=" ".join("""
            Path to sample diectories containing the crab base directories
            (see description).
        """.split()),
        metavar="PATH/TO/SAMPLE/DIRECTORIES",
        type=str,
        nargs="+",
    )

    args = parser.parse_args()
    if args.suffix == None:
        args.suffix = ["", "recovery_1", "recovery_2"]

    if args.status_files == None:
        args.status_files = ["status_0", "status_1", "status"]
    return args

if __name__ == '__main__':
    args = parse_arguments()
    main(**vars(args))