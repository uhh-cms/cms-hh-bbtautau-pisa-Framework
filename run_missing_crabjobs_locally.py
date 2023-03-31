import os
import sys
import json

from tqdm import tqdm
from datetime import datetime
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from subprocess import call

thisdir = os.path.realpath(os.path.dirname(__file__))

if not thisdir in sys.path:
    sys.path.append(thisdir)

from wlcg_dbs_interface import WLCGInterface

interface = WLCGInterface()

verbosity = 0

def run_custom_nano_command(
    input_file: str,
    nevents: int=100,
    era: str="Run2_2017"
):
    # first build the cms PSet
    cmd = [f'python3 {thisdir}/RunKit/nanoProdWrapper.py']
    cmd += [f'customise=Framework/NanoProd/customiseNano.customise'] 
    cmd += [f'skimCfg={thisdir}/config/skim.yaml']
    cmd += [f'maxEvents={nevents}'] 
    cmd += [f'sampleType=mc'] 
    cmd += [f'storeFailed=True'] 
    cmd += [f'era={era}'] 
    cmd += [f'inputFiles=file:{input_file}']
    print(f"creating PSet.py for input '{input_file}'")
    call([" ".join(cmd)], shell=True)
    cmd=f'{thisdir}/RunKit/nanoProdCrabJob.sh'

    print("executing job")
    call([cmd], shell=True)

def run_job(
    lfn: str,
    tmp_dir: str,
    wlcg_path:str,
    output_name:str,
    fail_on_exception: bool=False,
    **kwargs,
):
    # cd into tmp dir. If it doesn't exist, create it
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    cwd = os.getcwd()
    os.chdir(tmp_dir)

    local_lfn_name = os.path.basename(lfn)
    # copy lfn locally
    interface.get_remote_file(
        filepath=lfn,
        target=local_lfn_name
    )
    if not os.path.exists(local_lfn_name):
        msg = f"Unable to load '{lfn}'"
        if fail_on_exception:
            raise ValueError(msg)
        else:
            print(msg)
            os.chdir(cwd)
            return

    # run the job with this LFN
    run_custom_nano_command(
        input_file=local_lfn_name,

    )
    # copy the output to the (remote) WLCG site
    # from IPython import embed; embed()
    interface.move_file_to_remote(
        local_file=os.path.abspath("nano.root"),
        target_file=f"{wlcg_path}/{output_name}",
        route_url=None
    )

    # change back to the original directory after we're done
    os.chdir(cwd)


def build_wlcg_path(
    wlcg_prefix: str,
    wlcg_dir: str,
    sample_name: str,
    crab_dirname: str,
    time_stamp: str,
    job_output: str,
    wlcg_template: str=os.path.join("{wlcg_prefix}{wlcg_dir}",
        "{sample_name}",
        "{crab_dirname}",
        "{time_stamp}",
        "{job_output}"
    ),
):
    return wlcg_template.format(
        wlcg_prefix=wlcg_prefix,
        wlcg_dir=wlcg_dir,
        sample_name=sample_name,
        crab_dirname=crab_dirname,
        time_stamp=time_stamp,
        job_output=job_output,
    )

def main(
    *args, 
    wlcg_prefix: str,
    wlcg_dir: str,
    missing_files_json: str,
    sample_config: str,
    veto_dirs: list[str]=None,
    tmp_dir: str="./tmp",
    remote_dir: str="manual_jobs",
    **kwargs,
):
    if not veto_dirs:
        veto_dirs=list()
    # load dictionary with missing lfns
    missing_lfn_dict = dict()
    with open(missing_files_json) as f:
        missing_lfn_dict = json.load(f)
    
    # filter out samples that are accounted for in the list of veto directories
    missing_samples = list(filter(
        lambda x: not any(dir.strip("/").endswith(x) for dir in veto_dirs),
        missing_lfn_dict
    ))
    # from IPython import embed; embed()
    pbar_samples = tqdm(missing_samples)
    # create a time stamp that mimics the crab format
    timestamp = '{:%y%m%d_%H%M%S}'.format(datetime.now())
    local_job_summary = dict()
    for sample in pbar_samples:
        pbar_samples.set_description(f"Run missing jobs for sample '{sample}'")

        # load das key for this sample
        das_key = interface.load_das_key(
            sample_name=sample,
            sample_config=sample_config
        )

        # get campaign name
        sample_campaign = interface.get_campaign_name(das_key=das_key)
        
        # loop through the list of missing lfns
        pbar_missing_lfns = tqdm(missing_lfn_dict[sample]["missing_lfns"])
        missing_lfns = list()
        for i, lfn in enumerate(pbar_missing_lfns):
            lfn_shortname = "/".join(lfn.split("/")[-3:])
            pbar_missing_lfns.set_description(f"Running LFN {lfn_shortname}")

            blocknumber = int(i/10000)
            wlcg_path = build_wlcg_path(
                wlcg_prefix=wlcg_prefix,
                wlcg_dir=wlcg_dir,
                sample_name=sample_campaign,
                crab_dirname=remote_dir,
                time_stamp=timestamp,
                job_output=f"{blocknumber:04d}"
            )
            run_job(
                lfn=lfn,
                tmp_dir=os.path.join(tmp_dir, sample),
                wlcg_path=wlcg_path,
                output_name=f"nano_{i}.root",
            )
            missing_lfns.append(lfn)
        local_job_summary[sample] = {
            "timestamp": timestamp,
            "remote_dir": remote_dir,
            "lfns": missing_lfns,
        }
    with open("local_job_summary.json", "w") as f:
        json.dump(local_job_summary, f, indent=4)

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
        "-r", "--remote-dir",
        help=" ".join(
            """
                safe the output of the jobs you run locally in this
                directory in the (remote) WLCG site. The path on the remote
                site is formed as 
                '{wlcg-dir}/{sample_campaign}/{remote-dir}/{time-stamp}'.
                Defaults to 'manual_jobs'
            """.split()
        ),
        type=str,
        default="manual_jobs",
        dest="remote_dir",
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

    parser.add_argument(
        "veto_dirs",
        help=" ".join(
            """
                veto samples in the missing files json that are accounted
                for in this list. 
            """.split()
        ),
        metavar="path/to/crab_dirs",
        nargs="+"
    )

    args = parser.parse_args()
    # from IPython import embed; embed()
    
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