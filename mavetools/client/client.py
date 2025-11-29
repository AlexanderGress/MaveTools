import json
import logging
import requests
import sys
import os

from mavetools.models.scoreset import ScoreSet
from mavetools.models.ml_tools import MlExperiment


FUNCTION_TYPES = set(['Activity', 'Binding', 'Expression', 'OrganismalFitness', 'Stability'])

def extract_function_type(scoreset):
    function_type = None
    sd = scoreset['shortDescription']
    methodText = scoreset['methodText']
    abstract = scoreset['abstractText']
    if abstract is None or abstract == '':
        try:
            abstract = scoreset['primaryPublicationIdentifiers'][0]['abstract']
            if abstract is None:
                abstract = ''
        except IndexError:
            abstract = ''
        second_abstract = ''
    else:
        try:
            second_abstract = scoreset['primaryPublicationIdentifiers'][0]['abstract']
            if second_abstract is None:
                second_abstract = ''
        except IndexError:
            second_abstract = ''

    prio_binding_keywords = ['transcription ability', 'protein-protein interaction']
    prio_fitness_keywords = []
    prio_act_keywords = ['readout of activity', 'fluorescence of Aequorea victoria GFP', 'loss-of-function', 'Functional scores', 'functional scores']
    prio_exp_keywords = ['variant abundance', 'protein abundance']
    prio_stabi_keywords = []

    binding_keywords = ['binding', 'Binding', 'Y2H assay']
    fitness_keywords = ['fitness', 'Fitness', 'pathogenicity', 'saturation prime editing', 'Growth', 'growth', 'toxicity assay', 'TileSeq', 'viral replication', 'complementation', 'based on survival']
    act_keywords = ['activity', 'Activity']
    exp_keywords = ['expression', 'Expression', ]
    stabi_keywords = ['stability', 'Stability', 'folding free energy', 'dG value', 'Protein folding']

    keyword_tuples = [
        ('OrganismalFitness', fitness_keywords, prio_fitness_keywords),
        ('Binding', binding_keywords, prio_binding_keywords),
        ('Activity', act_keywords, prio_act_keywords),
        ('Expression', exp_keywords, prio_exp_keywords),
        ('Stability', stabi_keywords, prio_stabi_keywords)
    ]

    for ft, keywords, prio_keywords in keyword_tuples:
        for keyword in prio_keywords:
            if sd.count(keyword) > 0:
                function_type = ft
            elif methodText.count(keyword) > 0:
                function_type = ft
            if function_type is not None:
                break
        if function_type is not None:
            break

    if function_type is None:
        for ft, keywords, prio_keywords in keyword_tuples:
            for keyword in prio_keywords:
                if abstract.count(keyword) > 0:
                    function_type = ft
                    break
                elif second_abstract.count(keyword) > 0:
                    function_type = ft
                    break
            if function_type is not None:
                break

    if function_type is None:
        for ft, keywords, prio_keywords in keyword_tuples:
            for keyword in keywords:
                if sd.count(keyword) > 0:
                    function_type = ft
                elif methodText.count(keyword) > 0:
                    function_type = ft
                if function_type is not None:
                    break
            if function_type is not None:
                break

    if function_type is None:
        for ft, keywords, prio_keywords in keyword_tuples:
            for keyword in keywords:
                if abstract.count(keyword) > 0:
                    function_type = ft
                    break
                elif second_abstract.count(keyword) > 0:
                    function_type = ft
                    break
            if function_type is not None:
                break

    
    if function_type is not None:
        return function_type
    else:
        return 'Activity'

class ClientTemplate:
    """
    Parent class for client classes to inheir.
    """

    def parse_json_scoreset_list(
        self,
        scoreset_list,
        keywords=None,
        organisms=None,
        retrieve_json_only=False,
        experiment_types=None,
        verbose=False,
    ):
        """
        Parses a list of scoreset metadatas in json format.
        Creates classes and datastructures that are required to use ML tools features.

        Parameters
        ----------

        scoreset_list
            A list of scoreset metadata in json format.

        keywords
            List of keywords. If not None, filters all scoresets that keywords do not contain any of the given keywords.

        organisms
            List of organisms. If not None, filters all scoresets that organism is not any of the given organisms.

        retrieve_json_only
            When True, the function does not create the ML datastructures, but applies all given filters and creatures a ordered experiment-scoreset datastructure of the json objects.

        experiment_types
            List of experiment types. If not None, filters all scoresets that experiment type is not any of the given experiment types.

        Returns
        -------

        experiment_dict
            A dictionary mapping experiment urns to their corresponding MLExperiment objects.
        """

        if verbose:
            print(f"Call of parse_json_scoreset_list: {experiment_types=}")

        if keywords is not None:
            keywords = set(keywords)

        if organisms is not None:
            organisms = set(organisms)

        if experiment_types is not None:
            experiment_types = set(experiment_types)

        experiment_dict = {}
        filter_0 = 0
        filter_1 = 0
        filter_2 = 0
        filter_3 = 0
        filter_4 = 0
        for scoreset in scoreset_list:
            if keywords is not None:
                keyword_match = False
                for keyword in scoreset["keywords"]:
                    if keyword["text"] in keywords:
                        keyword_match = True
                if not keyword_match:
                    filter_0 += 1
                    continue

            if organisms is not None:
                if scoreset["target"]["reference_maps"][0]["genome"]["organism_name"] not in organisms:
                    filter_1 += 1
                    continue

            if len(scoreset["targetGenes"]) == 0:
                filter_2 += 1
                continue

            if experiment_types is not None:
                if scoreset["targetGenes"][0]["category"] not in experiment_types:
                    # print(scoreset['targetGenes'][0]['category'])
                    filter_3 += 1
                    continue

            urn = scoreset["urn"]

            if retrieve_json_only:
                experiment_dict[urn] = scoreset
                filter_4 += 1
                continue

            ft = extract_function_type(scoreset)

            scoreset_obj = ScoreSet.deserialize(scoreset)

            experiment_urn = scoreset_obj.urn

            if experiment_urn not in experiment_dict:
                experiment_dict[experiment_urn] = MlExperiment(
                    experiment_urn, {}, scoreset_obj, urn=experiment_urn, function_type = ft
                )

            experiment_dict[experiment_urn].scoreset_dict[urn] = scoreset_obj

        if verbose:
            print(
                f"Filtered: {filter_0=} {filter_1=} {filter_2=} {filter_3=} {filter_4=}"
            )

        return experiment_dict


class LocalClient(ClientTemplate):
    """
    A client class that imitates the original client class to use a local clone of the MaveDB.
    """

    def __init__(self, local_instance_path):
        """
        Initializes the client instance.

        Parameters
        ----------

        local_instance_path
            path the locally stored MaveDB.
        """

        self.local_instance_path = local_instance_path
        self.meta_data_folder = f"{local_instance_path}/main.json"
        self.main_meta_data = self.load_meta_data(self.meta_data_folder)
        self.scoreset_data_folder = f"{local_instance_path}/csv/"

    def get_meta_file_path(self, urn):
        """
        Getter for filepath of the stored meta data file in the locally cloned MaveDB.

        Parameters
        ----------

        urn
            MaveDB urn identifier of the scoreset.

        Returns
        -------

        path
            Path to the metadata json-formatted file.
        """
        return f"{self.meta_data_folder}/{urn}.json"

    def load_meta_data(self, filepath):
        """
        Wrapper function for loading a json-formatted metadata file.

        Parameters
        ----------

        filepath
            Path to a json-formatted metadata file.

        Returns
        -------

        meta_data
            json object of the metadata.
        """

        f = open(filepath, "r")
        meta_data = json.load(f)
        f.close()
        return meta_data

    def get_meta_data(self, urn):
        """
        Getter for metadata.

        Parameters
        ----------

        urn
            MaveDB urn identifier of a scoreset.

        Returns
        -------

        meta_data
            json object of the metadata.
        """
        return self.load_meta_data(self.get_meta_file_path(urn))

    def search_database(
        self,
        keywords=None,
        organisms=None,
        experiment_types=["protein_coding"],
        verbose=False,
    ):
        """
        Searches all scoresets in MaveDB and applies some filters.
        Parameters
        ----------

        keywords
            List of keywords. If not None, filters all scoresets that keywords do not contain any of the given keywords.

        organisms
            List of organisms. If not None, filters all scoresets that organism is not any of the given organisms.

        experiment_types
            List of experiment types. If not None, filters all scoresets that experiment type is not any of the given experiment types.

        Returns
        -------

        experiment_dict
            A dictionary mapping experiment urns to their corresponding MLExperiment objects.
        """

        experiment_sets = self.main_meta_data["experimentSets"]
        scoreset_list = []
        for experiment_set in experiment_sets:
            for experiment in experiment_set["experiments"]:
                for scoreSet in experiment["scoreSets"]:
                    scoreset_list.append(scoreSet)

        if verbose:
            print(f"Searching MaveDB: {len(experiment_sets)=} {len(scoreset_list)=}")

        experiment_dict = self.parse_json_scoreset_list(
            scoreset_list,
            keywords=keywords,
            organisms=organisms,
            experiment_types=experiment_types,
            verbose=verbose,
        )

        if verbose:
            print(f"{len(experiment_dict)=}")

        return experiment_dict

    def get_experiment_dict(self, urns):
        """
        Generates a experiment_dict containing MLExperiment objects for a list of given urns.

        Parameters
        ----------

        urns
            A list of MaveDB urn identifiers.

        Returns
        -------

        experiment_dict
            A dictionary mapping experiment urns to their corresponding MLExperiment objects.
        """

        scoreset_list = []
        for urn in urns:
            try:
                scoreset_list.append(self.get_meta_data(urn))
            except:
                n = 1
                while True:
                    try:
                        scoreset_urn = f"{urn}-{n}"
                        scoreset_list.append(self.get_meta_data(scoreset_urn))
                    except:
                        break
                    n += 1
        return self.parse_json_scoreset_list(scoreset_list)

    def retrieve_score_table(self, urn):
        """
        Retrieves the score table for an urn.

        Parameters
        ----------

        urn
            MaveDB urn identifier of a scoreset.

        Returns
        -------

        text
            Scoreset table as a string.
        """

        fixed_urn = urn.replace(":", "-")

        score_table_file = f"{self.scoreset_data_folder}/{fixed_urn}.scores.csv"
        f = open(score_table_file, "r")
        text = f.read()
        f.close()
        return text


class Client(ClientTemplate):
    def __init__(self, base_url="https://www.mavedb.org/api/", auth_token=""):
        """
        Instantiates the Client object and sets the values for base_url and
        auth_token

        Parameters
        ----------
        base_url: the url in which the api endpoint exists
            default: 'http://127.0.0.1:8000/api/'
        auth_token: authorizes POST requests via the API and MaveDB
            default: ''
        """
        self.base_url = base_url
        if auth_token:
            self.auth_token = auth_token

    class AuthTokenMissingException(Exception):
        pass

    def clone(self, local_instance_path):
        """
        Downloads the whole MaveDB and creates a local clone.

        Parameters
        ----------

        local_instance_path
            Path to where the clone should be stored.
        """

        if not os.path.exists(local_instance_path):
            os.mkdir(local_instance_path)
        meta_data_folder = f"{local_instance_path}/meta_data/"
        if not os.path.exists(meta_data_folder):
            os.mkdir(meta_data_folder)
        scoreset_data_folder = f"{local_instance_path}/scoreset_data/"
        if not os.path.exists(scoreset_data_folder):
            os.mkdir(scoreset_data_folder)

        entry_dict = self.search_database(retrieve_json_only=True)
        for urn in entry_dict:
            meta_file = f"{meta_data_folder}/{urn}.json"

            f = open(meta_file, "w")
            json.dump(entry_dict[urn], f)
            f.close()

            score_table_file = f"{scoreset_data_folder}/{urn}.csv"

            f = open(score_table_file, "w")
            f.write(self.retrieve_score_table(urn))
            f.close()

    def search_database(
        self,
        keywords=None,
        organisms=None,
        retrieve_json_only=False,
        experiment_types=["protein_coding"],
    ):
        """
        Searches all scoresets in MaveDB and applies some filters.
        Parameters
        ----------

        keywords
            List of keywords. If not None, filters all scoresets that keywords do not contain any of the given keywords.

        organisms
            List of organisms. If not None, filters all scoresets that organism is not any of the given organisms.

        experiment_types
            List of experiment types. If not None, filters all scoresets that experiment type is not any of the given experiment types.

        Returns
        -------

        experiment_dict
            A dictionary mapping experiment urns to their corresponding MLExperiment objects.
        """

        search_page_url = f"{self.base_url}/scoresets"
        try:
            r = requests.get(search_page_url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.json())
            raise SystemExit(e)

        scoreset_list = r.json()

        experiment_dict = self.parse_json_scoreset_list(
            scoreset_list,
            keywords=keywords,
            organisms=organisms,
            retrieve_json_only=retrieve_json_only,
        )

        return experiment_dict

    def retrieve_score_table(self, urn):
        """
        Retrieves the score table for an urn.

        Parameters
        ----------

        urn
            MaveDB urn identifier of a scoreset.

        Returns
        -------

        text
            Scoreset table as a string.
        """

        base_parent = self.base_url.replace("api/", "")
        score_table_url = f"{base_parent}scoreset/{urn}/scores/"
        try:
            r = requests.get(score_table_url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.json())
            raise SystemExit(e)
        return r.text

    def get_model_instance(self, model_class, instance_id):
        """
        Using a GET, hit an API endpoint to get info on a particular instance
        of a model class such as a ScoreSet.
        This will perform the HTTP GET request and then let the class itself
        parse the JSON data.

        Parameters
        ----------
        model_class : ModelClass
            The model class we want to which we want to cast the response.
            (e.g., Experiment or Scoreset)
        instance_id : str
            The id of the object we are retrieving.

        Returns
        -------
        model_instance
            An instance of the passed class.

        Raises
        ------
        ValueError
            If any mandatory fields are missing.
        """
        model_url = f"{self.base_url}{model_class.api_url()}"
        instance_url = f"{model_url}{instance_id}/"
        try:
            r = requests.get(instance_url)
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.json())
            raise SystemExit(e)
        return model_class.deserialize(r.json())

    def post_model_instance(self, model_instance):
        """
        Using a POST, hit an API endpoint to post a resource.
        Performs HTTP POST request.

        Parameters
        ----------
        model_instance
            instance of model that will be POSTed

        Returns
        -------
        requests.model.Response
            The HTTP response object from the request, which contains the URN
            of the newly-created model in the `Response.text` field.

        Raises
        ------
        AuthTokenMissingException
            If the auth_token is missing
        """

        # save object type of model_instance
        model_class = type(model_instance)
        model_url = f"{self.base_url}{model_class.api_url()}/"
        payload, files = model_instance.post_payload()

        # check for existance of self.auth_token, raise error if does not exist
        if not self.auth_token:
            error_message = "Need to include an auth token for POST requests!"
            logging.error(error_message)
            raise AuthTokenMissingException(error_message)

        try:  # to post data
            r = requests.post(
                model_url,
                data={"request": json.dumps(payload)},
                files=files,
                headers={"Authorization": (self.auth_token)},
            )
            r.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logging.error(r.text)
            sys.exit(1)

        # No errors or exceptions at this point, log successful upload
        logging.info(f"Successfully uploaded {model_instance}!")

        # return the HTTP response
        return r
