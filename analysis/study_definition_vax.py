from cohortextractor import (
    StudyDefinition, 
    patients, 
    filter_codes_by_category
)

# Import codelists.py script
from codelists import *

import pandas as pd

# import the variables for deriving JCVI groups
from grouping_variables import (
    jcvi_variables, 
    study_parameters
)

# define variables explicitly from study_parameters
max_comparisons=study_parameters["K"]
start_date=study_parameters["start_date"] # start of phase 1
end_date=study_parameters["end_date"] # latest date of data

# regions
regions = pd.read_csv(
    filepath_or_buffer='./analysis/lib/regions.csv',
    dtype=str
)
ratio_regions = { regions['region'][i] : float(regions['ratio'][i]) for i in regions.index }

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

# Will's function for extracting vaccination dates
def vaccination_date_X(name, index_date, n, product_name_matches=None, target_disease_matches=None):
  # vaccination date, given product_name
  def var_signature(
    name,
    on_or_after,
    product_name_matches,
    target_disease_matches
  ):
    return {
      name: patients.with_tpp_vaccination_record(
        product_name_matches=product_name_matches,
        target_disease_matches=target_disease_matches,
        on_or_after=on_or_after,
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
      ),
    }
    
  variables = var_signature(f"{name}_1_date", index_date, product_name_matches, target_disease_matches)
  for i in range(2, n+1):
    variables.update(var_signature(
      f"{name}_{i}_date", 
      f"{name}_{i-1}_date + 1 days",
      # pick up subsequent vaccines occurring one day or later -- people with unrealistic dosing intervals are later excluded
      product_name_matches,
      target_disease_matches
    ))
  return variables

study=StudyDefinition(

    default_expectations={
        "date": {"earliest": start_date, "latest": end_date},
        "rate": "uniform",
        "incidence": 0.8,
    },

    **jcvi_variables,   

    population=patients.satisfying(
        """
        has_follow_up = 1
        AND
        NOT died_before
        """,
        has_follow_up=patients.registered_with_one_practice_between(
            start_date="elig_date - 1 year",
            end_date="elig_date",
            return_expectations={"incidence": 0.90},
        ),
        died_before=patients.died_from_any_cause(
            on_or_before="elig_date + 42 days",
            returning="binary_flag",
        ),
    ),

    # Healthcare worker flag on vaccine record
    hscworker=patients.with_healthcare_worker_flag_on_covid_vaccine_record(
        returning="binary_flag",
        return_expectations={"incidence": 0.01},
        ),

    ######################
    ### COVID VACCINES ###
    ######################

    # pfizer
  **vaccination_date_X(
    name = "covid_vax_pfizer",
    # use 1900 to capture all possible recorded covid vaccinations, including date errors
    # any vaccines occurring before national rollout are later excluded
    index_date = "1900-01-01", 
    n = 3,
    product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)"
  ),
  
  # az
  **vaccination_date_X(
    name = "covid_vax_az",
    index_date = "1900-01-01",
    n = 3,
    product_name_matches="COVID-19 Vaccine Vaxzevria 0.5ml inj multidose vials (AstraZeneca)"
  ),
  
  # moderna
  **vaccination_date_X(
    name = "covid_vax_moderna",
    index_date = "1900-01-01",
    n = 3,
    product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)"
  ),
  
  # any covid vaccine
  # use n=1 here as this is just used to exclude individuals from the unvaccinated arm
  **vaccination_date_X(
    name = "covid_vax_disease",
    index_date = "1900-01-01",
    n = 1,
    target_disease_matches="SARS-2 CORONAVIRUS"
  ),

    #############################
    ### DEMOGRAPHIC VARIABLES ###
    #############################

    # ETHNICITY IN 6 CATEGORIES
    # ethnicity
    ethnicity_6=patients.with_these_clinical_events(
        eth2001_primis,
        returning="category",
        find_last_match_in_period=True,
        on_or_before="elig_date + 42 days",
        return_expectations={
            "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
            "incidence": 0.75,
        },
    ),

    ethnicity_6_sus=patients.with_ethnicity_from_sus(
        returning="group_6",  
        use_most_frequent_code=True,
        return_expectations={
            "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
            "incidence": 0.8,
        },
    ),

    # Patients in long-stay nursing and residential care
    # any time before end_K_date (return earliest date)
    longres_date=patients.with_these_clinical_events(
        longres_primis,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_after=start_date,
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),

    # region, IMD and medically housebound defined on or before elig_date + 42 days,
    # as this is the time-point at which they are used in eligibility criteria
    # region - NHS England 9 regions
    region=patients.registered_practice_as_of(
        "elig_date + 42 days",
        returning="nuts1_region_name",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": ratio_regions,
            },
            "incidence": 0.99
        },
    ),

    # IMD
    imd=patients.address_as_of(
                    "elig_date + 42 days",
                    returning="index_of_multiple_deprivation",
                    round_to_nearest=100,
                    return_expectations={
                        "category": {"ratios": {c: 1/320 for c in range(100,32100,100)}}
                        }
                    ),

    # medically housebound
    housebound = patients.satisfying(
        """
        housebound_date
        AND NOT no_longer_housebound
        AND NOT moved_into_care_home
        """,
        
    housebound_date=patients.with_these_clinical_events( 
        housebound, 
        on_or_before="elig_date + 42 days",
        find_last_match_in_period = True,
        returning="date",
        date_format="YYYY-MM-DD",
    ),   
    no_longer_housebound=patients.with_these_clinical_events( 
        no_longer_housebound, 
        on_or_after="housebound_date",
    ),
    moved_into_care_home=patients.with_these_clinical_events(
        longres_primis,
        on_or_after="housebound_date",
        ),
    ),
    
    # End of life care defined on or before elig_date +84 days, as this will be defined
    # at the start of the SVP, and this is the latest date at which the SVP can start.
    endoflife_date=patients.with_these_clinical_events(
        eol_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 84 days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    # midazolam
    midazolam_date=patients.with_these_medications(
        midazolam_codes,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before="elig_date + 84 days",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "exponential_increase",
            "incidence": 0.01
        },
    ),
    
    ##############
    ### EVENTS ###
    ##############
    
    # any death
    death_date=patients.died_from_any_cause(
        returning="date_of_death",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date},
            "rate": "uniform",
            "incidence": 0.02
        },
    ),
    # De-registration
    dereg_date=patients.date_deregistered_from_all_supported_practices(
        on_or_after="elig_date",
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": start_date, "latest": end_date,},
            "incidence": 0.001
        }
    ),

)
