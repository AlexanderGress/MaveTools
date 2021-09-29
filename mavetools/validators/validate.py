from mavetools.validators import dataset_validators

def validate_all():
    '''
    By calling other helper functions, this function runs all of the validation code
    '''
    return None

def validate_dataset(scorefile=None, countfile=None, scorejson=None):
    '''
    This function calls all of the validation functions within
    mavetools/mavetools/validators/dataset_validation.py
    Returns
    -------

    '''

    if scorefile is not None:
        # this one returns header
        scoreheader = dataset_validators.read_header_from_io(file=scorefile)

        # if the header was returned, do these ones
        dataset_validators.validate_has_hgvs_in_header(header=scoreheader)
        dataset_validators.validate_at_least_one_additional_column(header=scoreheader)
        dataset_validators.validate_header_contains_no_null_columns(header=scoreheader)

        dataset_validators.validate_scoreset_score_data_input(file=scorefile)

        dataset_validators.validate_scoreset_json(dict_=scorejson)

    if countfile is not None:
        countheader = dataset_validators.read_header_from_io(file=countfile)

        # if the header was returned, do these ones
        dataset_validators.validate_has_hgvs_in_header(header=countheader)
        dataset_validators.validate_at_least_one_additional_column(header=countheader)
        dataset_validators.validate_header_contains_no_null_columns(header=countheader)

        dataset_validators.validate_scoreset_count_data_input(file=countfile)

    if scorefile is not None and countfile is not None:
        dataset_validators.validate_datasets_define_same_variants(scores=scorefile, counts=countfile)




