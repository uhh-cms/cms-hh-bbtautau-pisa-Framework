import os
import sys
import yaml
import json

try:
    import gfal2
except ImportError as e:
    print("WARNING: could not import gfal2")
    print(e)
    print("gfal will be disabled!")
    gfal2 = None



from subprocess import PIPE, Popen
from itertools import chain
from typing import Any
from tqdm import tqdm

class WLCGInterface(object):
    def __init__(self,
        # wlcg_path: str or None=None,
        # route_url: str or None=None,
        verbosity: int=0,
    ):
        # self.wlcg_path = wlcg_path
        # self.route_url = route_url
        self.__verbosity = verbosity
        # setup gfal context
        try:
            # create gfal context
            if not gfal2:
                raise NotImplementedError("Cannot load remote file without gfal2 module!")

            self.gfal_context = gfal2.creat_context()
        except NotImplementedError as e:
            print(e)
            self.gfal_context = None
        self.dbs_api = self.setup_dbs_api()


    @property
    def verbosity(self):
        return self.__verbosity
    

    @verbosity.setter
    def verbosity(self, val: int):
        self.__verbosity = val

    
    def setup_dbs_api(
        self,
        cms_dbs_url: str="https://cmsweb.cern.ch/dbs/prod/global/DBSReader",
    ):
        # setup dbs search
        try:
            from dbs.apis.dbsClient import DbsApi
            
            return DbsApi(url=cms_dbs_url)
        except:
            print("WARNING: Could not find dbs3 module. Did you install it with")
            print("python3 -m pip install --user dbs3-client")
            print("?")
            print("Will use dasgoclient as fallback instead")
            return None

    def get_remote_file(
        self,
        filepath:str ,
        target: str,
        route_url: str="root://cms-xrd-global.cern.ch",
    ):  
        url=f"{route_url}//{filepath}"
        target = os.path.abspath(target)

        # copy remote file
        self.gfal_context.filecopy(
            self.gfal_context.transfer_parameters(),
            url,
            f"file://{target}"
        )


    def load_remote_output(
        self,
        wlcg_path: str,
    ) -> list[str]:
        """Function to load file paths from a remote WLCG target *wlcg_path*.
        First, the function checks for the gfal2 module. If gfal is loaded correctly,
        the list of files from the remote directly *wlcg_path* is loaded.
        If any of these steps fail, an empty list is returned

        Args:
            wlcg_path (str):    Path to the WLCG remote target, which consists of the
                                WLCG prefix and the actual directory on the remote
                                site (constructed from global wlcg_template)

        Returns:
            list[str]:  List of files in remote target *wlcg_path*. Defaults to
                        empty list.
        """
        try:
            if self.gfal_context:
                # load list of files
                filelist = self.gfal_context.listdir(wlcg_path)
                return [os.path.join(wlcg_path, x) for x in filelist]
            else:
                if self.verbosity >= 1:
                    print(f"unable to load files from {wlcg_path}, skipping")
                    from IPython import embed; embed()
        except Exception as e:
            print(f"unable to load files from {wlcg_path}, skipping")
        return []

# def load_events()

    def compare_events(
        self,
        relevant_ids,
        job_outputs,
        input_map,
        event_lookup,
    ):
        pbar_ids = tqdm(relevant_ids)
        for id in pbar_ids:
            pbar_ids.set_description(f"Comparing events for job {id}")
            relevant_job_outputs = set()
            relevant_job_outputs = set(filter(
                    lambda x: x.endswith(f"nano_{id}.root"), 
                    job_outputs
                ))
            all_events = sum([event_lookup.get(x, 0) for x in input_map[id]])


    def check_job_outputs(
        self,
        collector_set: set[str],
        input_map: dict[str, list[str]],
        job_details: dict[str, dict],
        state: str="failed",
        job_outputs: set or None=None,
        event_lookup: dict or None=None,
        verbosity: int=0,
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
        If the set *job_outputs* is neither None nor empty, the file paths are
        matched to the job ids with the current state. Only IDs with output files
        are considered further.

        Args:
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
            job_outputs (set, optional):    if a set of output files is given,
                                            only job ids with output files are
                                            considered as relevant. Defaults to None

        Raises:
            ValueError: If a lfn is already marked as done but is associated with
                        a done job again, the ValueError is raised.
        """
        relevant_ids = set(filter(
            lambda x: job_details[x]["State"] == state,
            job_details
        ))

        # if there are paths to the job outputs available, only select ids that
        # actually have an output
        if isinstance(job_outputs, set) and not len(job_outputs) == 0:
            relevant_ids = set(filter(
                lambda x: any(path.endswith(f"nano_{x}.root") for path in job_outputs), 
                relevant_ids
            ))

        # for state "failed", collect output files that should not be there
        if state == "failed":
            collector_set.update(filter(
                lambda x: any(x.endswith(f"nano_{id}.root") for id in relevant_ids), 
                job_outputs
            ))
        # if state is finished, safe the done lfns (if the output of the job is also 
        # available)
        elif state == "finished":
            
            lfns = set()
            
            # first check if a lfn is already marked as done - this should not happen
            
            if event_lookup:
                self.compare_events()
            else:
                lfns = set(chain.from_iterable([input_map[x] for x in relevant_ids]))
            overlap = collector_set.intersection(lfns)
            if len(overlap) != 0:
                if verbosity == 0:
                    msg = " ".join(f"""
                        {len(overlap)} LFNs are already marked as 'done' but
                        have come up here again. This should not happen
                    """.split())
                    raise ValueError(msg)
                else:
                    overlap_string = "\n".join(overlap)
                    raise ValueError(f"""
                    The following lfns were already marked as done:
                    {overlap_string}

                    This should not happen!
                    """)
                
            collector_set.update(lfns)


    def load_das_key(
        self,
        sample_name: str,
        sample_config: str,
        verbosity: int=0
    ) -> str or None:
        """Small function to extract the DAS key for sample *sample_name* 
        from the *sample_config*. First, the *sample_config*
        is opened (has to be in yaml format!). Afterwards, the entry *sample_name*
        is extracted. This entry should be a dictionary itself, which should contain
        the key 'miniAOD' with the DAS key for this sample.

        Args:
            sample_name (str): Name of the sample as provided in the sample config
            sample_config (str):    path to the sample_config.yaml file containing 
                                    the information mentioned above.

        Returns:
            str or None: If successfull, this function returns the DAS key, else None
        """    
        das_key = None
        # open the sample config
        with open(sample_config) as f:
            sample_dict = yaml.load(f, yaml.Loader)

        # look up information for sample_name
        sample_info = sample_dict.get(sample_name, dict())
        # if there is no sample information, exit here
        if len(sample_info) == 0:
            if verbosity >= 1:
                print(f"WARNING: Unable to load information for sample '{sample_name}'")
            return das_key
        return sample_info.get("miniAOD", None)

    def get_campaign_name(self, das_key: str=None, verbosity: int=0) -> str:
        """small function to translate the sample name attributed by the 
        crabOverseer to the original MC campaign name. The original 
        campaign name is then extracted from the DAS key. If any of these steps
        fails, the fucntion returns an empty string

        Args:
            das_key (str):  DAS key in str format. Any other format will return
                            the default value of '""'.

        Returns:
            str: if successful, returns the original campaign name, else ""
        """    
        sample_campaign = ""
        # load information about the original miniAOD DAS key
        if not isinstance(das_key, str):
            if verbosity >= 1:
                msg=" ".join(f"""
                WARNING: Unable to load campaign name from das key of type
                '{type(das_key)}'
                """.split())
                print(msg)
            return sample_campaign

        # original campaign name is the first part of the DAS key
        sample_campaign = das_key.split("/")[1]
        
        return sample_campaign

    def load_valid_file_list(self, das_key: str) -> dict[str: Any]:
        # load the file list for this dataset
        file_list = self.dbs_api.listFiles(dataset=das_key, detail=1)
        # by default, this list contains _all_ files (also LFNs that are not
        # reachable) so filter out broken files
        file_list = list(filter(
            lambda x: x["is_file_valid"] == True,
            file_list
        ))
        return file_list
    
    def create_event_lookup(
        self,
        das_key: str
    ) -> dict[str, int]:
        if self.dbs_api:
            dbs_file_list = self.load_valid_file_list(das_key=das_key)
            return {
                x["logical_file_name"]: x["event_count"] for x in dbs_file_list
            }
        return dict()

    def get_dbs_lfns(self, das_key: str) -> set[str]:
        """Small function to load complete list of valid LFNs for dataset with
        DAS key *das_key*. Only files where the flag 'is_file_valid' is True
        are considered. Returns set of lfn paths if successful, else an empty set

        Args:
            das_key (str): key in CMS DBS service for the dataset of interest

        Returns:
            set[str]: Set of LFN paths
        """    
        # initialize output set as empty
        output_set = set()

        # if the api for the dbs interface was initialized sucessfully, we can
        # load the files
        if self.dbs_api:
            file_list = self.load_valid_file_list(das_key=das_key)
            output_set = set([x["logical_file_name"] for x in file_list])
        return output_set 

    def get_das_information(
        self,
        das_key: str,
        relevant_info: str="num_file",
        default: int=-1,
    ) -> int:
        allowed_modes="file_size num_event num_file".split()
        if not relevant_info in allowed_modes:
            raise ValueError(f"""Could not load information '{relevant_info}'
            because it's not part of the allowed modes: file_size, num_event, num_file
            """)
        output_value = default

        # execute DAS query for sample with *das_key*
        process = Popen(
            [f"dasgoclient --query '{das_key}' -json"], 
            shell=True, stdin=PIPE, stdout=PIPE
        )
        # load output of query
        output, stderr = process.communicate()
        # from IPython import embed; embed()
        # output is a string of form list(dict()) and can be parsed with
        # the json module
        try:
            das_infos = json.loads(output)
        except Exception as e:
            # something went wrong in the parsing, so just return the default
            return output_value
        # not all dicts have the same (relevant) information, so go look for the
        # correct entry in list. Relevant information for us is the total
        # number of LFNs 'nfiles'
        relevant_values = list(set(
            y.get(relevant_info) 
            # only consider entries in DAS info list with dataset information
            for x in das_infos if "dataset" in x.keys() 
            # only consider those elements in dataset info that also have 
            # the relevant information
            for y in x["dataset"] if relevant_info in y.keys()
        ))
        # if this set is has more than 1 or zero entries, something went wrong
        # when obtaining the relevant information, so return the default value

        # if the set has exactly one entry, we were able to extract the relevant
        # information, so return it accordingly
        if len(relevant_values) == 1:
            output_value = relevant_values[0]
        return output_value