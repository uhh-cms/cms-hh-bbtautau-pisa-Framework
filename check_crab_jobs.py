import os
import sys
import json

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from tqdm import tqdm
from itertools import chain
from collections.abc import Iterable
from typing import Any
from subprocess import call
import numpy as np

thisdir = os.path.realpath(os.path.dirname(__file__))

if not thisdir in sys.path:
    sys.path.append(thisdir)

from wlcg_dbs_interface import WLCGInterface
interface = WLCGInterface()
wlcg_template= os.path.join("{wlcg_prefix}{wlcg_dir}",
    "{sample_name}",
    "{crab_dirname}",
    "{time_stamp}",
)
verbosity=0


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
        print(f"Could not load input mapping from file '{path}'")
        return None

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
        msg = "Project dir in status file does not match current crab dir under scutiny!"
        raise ValueError(msg)

def check_crab_directory(
    sample_dir: str,
    sample_name: str,
    suffix: str,
    das_key: str,
    status_file: str,
    known_lfns: set[str],
    done_lfns: set[str],
    failed_job_outputs: set[str],
    pbar: Iterable,
    wlcg_dir: str,
    wlcg_prefix: str,
    xrd_prefix: str,
    time_stamps: list[str],
    job_input_file: str="job_input_files.json",
    event_lookup: dict[str, int] or None=None,
    event_comparison_container: list[dict[str, Any]] or None=None,
    **kwargs,
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
        sample_campaing (str):  Name of the official campaign for sample 
                                *sample_name*. Can be obtained with 
                                meth::`parse_sample_name` based on the
                                sample_name and the sample_config.yaml file
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
        if verbosity >= 1:
            print(f"Directory {crab_dir} does not exist, will stop looking here")
        return
    pbar.set_description(f"Checking directory {crab_dir}")
    
    # load the input file mapping of the form 'job_ids' -> list of lfns
    input_map = get_job_inputs(crab_dir=crab_dir, job_input_file=job_input_file)
    
    if not input_map:
        if verbosity >= 1:
            print(f"WARNING: could not load input map for directory {crab_dir}")
        return
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
        msg = ""
        if verbosity == 0:
            msg = f"""
            Encountered {len(unknown_lfns)} while processing dir '{crab_dir}'
            This is not expected. For more information, use higher level of
            verbosity.
            """
            # raise ValueError(msg)
            
        else:
            unknown_lfns_string = "\n".join(unknown_lfns)
            from IPython import embed; embed()
            # raise ValueError(f"""
            #     Following lfns are not known in '{crab_dir}':
            #     {unknown_lfns_string}

            #     This should not happen!
            #     """)
        known_lfns.update(flat_lfns)

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
    campaign_name = interface.get_campaign_name(
        das_key=das_key,
    )
    this_wlcg_template = wlcg_template.format(
        wlcg_prefix=wlcg_prefix,
        wlcg_dir=wlcg_dir,
        sample_name=campaign_name,
        crab_dirname=crab_dirname,
        time_stamp=time_stamp
    )

    # create complete path on remote WLCG system to output file
    # crab arranges the output files in blocks depending on the job id
    
    # get maximum ID to identify maximum block number later
    max_jobid = np.max([int(x) for x in job_details.keys()])
    # from IPython import embed; embed()
    # initialize set of outputs
    job_outputs = set()
    # get the maximum block number and iterate through the blocks
    pbar_blocks = tqdm(range(int(max_jobid/1000)+1))
    for i in pbar_blocks:
        pbar_blocks.set_description(f"Loading outputs for block {i:04d}")
        job_outputs.update(
            interface.load_remote_output(
                wlcg_path=os.path.join(this_wlcg_template, f"{i:04d}"),
            )
        )

    # load information about failed jobs
    interface.check_job_outputs(
        job_outputs=job_outputs,
        collector_set=failed_job_outputs,
        input_map=input_map,
        job_details=job_details,
        state="failed",
        wlcg_prefix=wlcg_prefix,
        xrd_prefix=xrd_prefix,
    )

    # load information about finished jobs
    ndone = len(done_lfns)
    interface.check_job_outputs(
        job_outputs=job_outputs,
        collector_set=done_lfns,
        input_map=input_map,
        job_details=job_details,
        state="finished",
        event_comparison_container=event_comparison_container,
        event_lookup=event_lookup,
    )

    if ndone< len(done_lfns):
        time_stamps.append(time_stamp)
    else:
        time_stamps.append([])

def post_processing(
    meta_infos: dict[str, Any],
    event_comparison: dict[str, list[dict[str, Any]]] or None=None
):
    # do some final sanity check: if we found missing lfns, print them here
    # so the user can do something
    build_meta_info_table(meta_infos=meta_infos)

    samples_with_missing_lfns = list(filter(
        lambda x: meta_infos[x]["missing"] != 0 or 
                   ( meta_infos[x]["das_total"] != meta_infos[x]["total"]
                        and meta_infos[x]["das_total"] != -1
                    ), 
        meta_infos
    ))
    if len(samples_with_missing_lfns) != 0:
        print("\n\n")
        print("Samples with missing LFNS:")
        build_meta_info_table(
            meta_infos={x: meta_infos[x] for x in samples_with_missing_lfns},
            outfilename="samples_with_missing_lfns.json"
        )

    samples_wo_missing_lfns = list(filter(
        lambda x: meta_infos[x]["missing"] == 0 or 
                   ( meta_infos[x]["das_total"] != meta_infos[x]["total"]
                        and meta_infos[x]["das_total"] != -1
                    ), 
        meta_infos
    ))
    if len(samples_wo_missing_lfns) != 0:
        print("\n\n")
        print("Samples w/o missing LFNS:")
        build_meta_info_table(
            meta_infos={x: meta_infos[x] for x in samples_wo_missing_lfns},
            outfilename="samples_wo_missing_lfns.json"
        )
    
    if event_comparison and len(event_comparison) > 0:
        # save event comparison
        with open("event_comparison.json", "w") as f:
            json.dump(event_comparison, f, indent=4)
    
    

def main(*args,
    sample_dirs=[],
    # wlcg_dir=None,
    # wlcg_prefix="",
    suffices=[],
    status_files=[],
    sample_config=None,
    dump_filelists=False,
    rm_failed=False,
    local_job_summary=None,
    **kwargs
):
    """main function. Load information provided by the ArgumentParser. Loops
    Thorugh the sample directories provided as *sample_dirs* and the *suffices*
    to check the individual crab base directories.
    Finally, check if any lfns are unaccounted for in the list of finished jobs.
    """  
    # load the information from the argument parser  
    meta_infos = dict()
    event_comparison = dict()

    local_job_summary_dict = dict()
    if local_job_summary:
        with open(local_job_summary) as f:
            local_job_summary_dict = json.load(f)
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
        das_key = interface.load_das_key(
            sample_name=sample_name, sample_config=sample_config,
        )

        # get full set of lfns for this sample

        # if verbosity is >= 2, we perform an event comparison, 
        # so create lookup map accordingly
        event_lookup = None
        # container for event comparisons
        sample_event_comparison = None
        sum_events = None
        if verbosity >= 1:
            event_lookup = interface.create_event_lookup(das_key=das_key)
            # the list of lfns is now the list of keys
            known_lfns = set(event_lookup.keys())
            sum_events = sum(event_lookup.values())
            if verbosity >= 2:
                sample_event_comparison = list()
            else:
                event_lookup = None
        else:
            # otherwise, there is no need to look up the events, so just 
            # create the set of lfns directly
            known_lfns=interface.get_dbs_lfns(das_key=das_key)   

        # if the dbs could not be contacted for some reason, use DAS
        # to load the total number of LFNS
        if len(known_lfns) > 0:
            n_total = len(known_lfns)
        else:
            # get total number of LFNs from DAS
            n_total = interface.get_das_information(
                das_key=das_key
            )

        # set up the sets to keep track of the lfns
        done_lfns=set()     # set of lfns processed by successful jobs
        
        # set of **outputs** from failed jobs, which shouldn't happen
        failed_job_outputs=set()    

        # set of relevant time stamps (needed for later merging of files)
        time_stamps = list()

        

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
                das_key=das_key,
                time_stamps=time_stamps,
                event_comparison_container=sample_event_comparison,
                event_lookup=event_lookup,
                **kwargs,
            )
        
        local_job_infos = local_job_summary_dict.get(sample_name)
        if local_job_infos:
            done_lfns.update(local_job_infos["lfns"])
            time_stamps.append(local_job_infos["timestamp"])

        # in the end, all LFNs should be accounted for
        unprocessed_lfns = known_lfns.symmetric_difference(done_lfns)

        sample_dict = dict()
        sample_dict["das_total"] = n_total
        sample_dict["total"] = len(known_lfns)
        if sum_events:
            sample_dict["sum_events"] = sum_events
        sample_dict["done"] = len(done_lfns)

        if rm_failed:
            # rm_bar = tqdm(failed_job_outputs)
            # for f in rm_bar:
                # rm_bar.set_description(f"Deleting file {f}")
            def chunks(lst, n):
                """Yield successive n-sized chunks from lst."""
                for i in range(0, len(lst), n):
                    yield lst[i:i + n]
            for c in list(chunks(list(failed_job_outputs), 400)):      
                cmd = f"gfal-rm {' '.join(c)}"
                call([cmd], shell=True)
            failed_job_outputs = set()

        if len(failed_job_outputs) > 0:
            sample_dict["outputs from failed jobs"] = len(failed_job_outputs)
        sample_dict["missing"] = len(unprocessed_lfns)
        sample_dict["time_stamps"] = time_stamps.copy()
        if dump_filelists:
            sample_dict["total_lfns"] = list(known_lfns.copy())
            sample_dict["done_lfns"] = list(done_lfns.copy())
            sample_dict["missing_lfns"] = list(unprocessed_lfns.copy())
            if len(failed_job_outputs) > 0:
                sample_dict["failed_outputs"] = list(failed_job_outputs.copy())
        meta_infos[sample_name] = sample_dict.copy()
        if sample_event_comparison and len(sample_event_comparison) > 0:
            event_comparison[sample_name] = sample_event_comparison.copy()
        elif len(unprocessed_lfns) == 0 and len(failed_job_outputs) > 0:
            # if the event comparison contains nothing, it might indicate
            # that we actually processed all events.
            # We could then consider to delete the job outputs that were 
            # generated from failed jobs
            pass
        
        if verbosity >= 3:
            if len(failed_job_outputs) != 0:
                print("WARNING: found job outputs that should not be there")
                print(f"Sample: {sample_dir}")
                for f in failed_job_outputs:
                    print(f)
            if len(unprocessed_lfns) != 0:
                print(f"WARNING: following LFNs for sample {sample_dir} were not processed!")
                for f in unprocessed_lfns:
                    print(f)
    
    post_processing(meta_infos=meta_infos, event_comparison=event_comparison)

def build_meta_info_table(
    meta_infos: dict,
    outfilename: str="crab_job_summary.json"
) -> None:
    """Helper function to convert collected information of crab jobs into
    human-readable table. The information is safed in *meta_infos*, which
    is a dictionary of the format
    {
        sample_name: {
            "das_total": TOTAL_NUMBER_OF_LFNS_IN_DAS
            "total": TOTAL_NUMBER_OF_LFNS,
            "done": NUMBER_OF_DONE_LFNS,
            "missing": NUMBER_OF_MISSING_LFNS,
            "outputs from failed jobs": NUMBER_OF_OUTPUTS_FROM_FAILED_JOBS
        }
    }

    Current supported formats: markdown (md)

    Args:
        meta_infos (dict): Dictionary containing above mentioned information
    """    

    # first build the header of the table
    headerparts = ["{: ^16}".format(x) 
                for x in ["Sample", "DAS Total #LFNs", "Total #LFNs", 
                            "Done #LFNs", "Missing #LFNs", 
                            "#Outputs from failed",
                        ]
            ]
    lines = ["| {} |".format(" | ".join(headerparts))]
    # markdown tables have a separation line with one '---' per column
    lines += ["| {} |".format(" | ".join(["---"]*len(headerparts)))]
    # loop through the samples that were processed
    for s in meta_infos:
        # load information for table
        das_tot = meta_infos[s]["das_total"]
        tot = meta_infos[s]["total"]
        done = meta_infos[s]["done"]
        missing = meta_infos[s]["missing"]
        failed = meta_infos[s].get("outputs from failed jobs", 0)
        # create line for table
        lines.append("| {} |".format(" | ".join(
                ["{: ^16}".format(x) 
                    for x in [s, das_tot, tot, done, missing, failed ]
                ]
            )
            )
        )
    # create final table
    table = "\n".join(lines)
    print(table)
    with open(outfilename, "w") as f:
        json.dump(meta_infos, f, indent=4)


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

    If you want to check the outputs of the jobs on the (remote) target site,
    make sure that gfal2 is available to python. On lxplus, you can do this
    by first using

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
        "--xrd-prefix",
        help=" ".join(
            """
                Prefix to contact the directory at the remote site via XROOTD.
                Defaults to prefix for T2_DESY 
                (root://dcache-cms-xrootd.desy.de:1094)
            """.split()
        ),
        type=str,
        default="root://dcache-cms-xrootd.desy.de:1094",
        dest="xrd_prefix"
    )
    parser.add_argument(
        "-s", "--suffices",
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
        "--dump-filelists",
        help=" ".join(
            """
            save the paths to all lfns, the done lfns, the missing lfns
            and to outputs on the WLCG remote site from failed jobs
            in the final summary .json file. Defaults to False
            """.split()
        ),
        default=False,
        action="store_true",
        dest="dump_filelists",
    )

    parser.add_argument(
        "--rm-failed",
        help=" ".join(
            """
            directly remove job outputs from failed jobs at the remote
            target location. WARNING: this is not reversible! Please make
            sure you know what you're doing!
            """.split()
        ),
        default=False,
        action="store_true",
        dest="rm_failed",
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

    parser.add_argument(
        "-v", "--verbosity",
        type=int,
        default=0,
        help=" ".join(
            """
                control the verbosity of the output. Currently implemented levels
                0:  just print number of files/outputs for each sample
                2:  actually check the event contents of the outputs w.r.t. the
                    corresponding LFN content
                3:  actually print the paths for the different files

            """.split()
        )
    )


    parser.add_argument("-l", "--local-job-summary",
        help=" ".join(
            """
            path to summary files about locally run jobs
            """.split()
        ),
        metavar="path/to/summary_for_local_jobs.json",
        # nargs="+",
        type=str,
        dest="local_job_summary"
    )

    args = parser.parse_args()
    if args.suffices == None:
        args.suffices = ["", "recovery_1", "recovery_2"]

    if args.status_files == None:
        args.status_files = ["status_0", "status_1", "status"]
    
    if not os.path.exists(args.sample_config):
        parser.error(f"file {args.sample_config} does not exist!")
    
    global verbosity
    verbosity = args.verbosity
    interface.verbosity = verbosity
    return args

if __name__ == '__main__':
    args = parse_arguments()
    main(**vars(args))