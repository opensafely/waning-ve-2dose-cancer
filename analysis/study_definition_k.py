from cohortextractor import (
    StudyDefinition, 
    patients
)

# Import codelists.py script
from codelists import *

# import json module
import json

#study_parameters
with open("./analysis/lib/study_parameters.json") as f:
  study_parameters = json.load(f)

# define variables explicitly from study_parameters
max_comparisons=study_parameters["K"]
start_date=study_parameters["start_date"] # start of phase 1
end_date=study_parameters["end_date"] # latest date of data

# set seed so that dummy data can be reproduced
import numpy as np
np.random.seed(study_parameters["seed"])

k="%placeholder_k%"

###
study=StudyDefinition(

    default_expectations={
        "date": {"earliest": start_date, "latest": end_date},
        "rate": "uniform",
        "incidence": 0.8,
    },  

    population=patients.all(),

    # start date of comparison period k
    start_k_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning=f'start_{k}_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),
    # end date of comparison period k (start +28 days for vax, +56 days for unvax)
    end_k_date=patients.with_value_from_file(
        f_path='output/data/data_eligible_e.csv', 
        returning=f'end_{k}_date', 
        returning_type='date',
        date_format='YYYY-MM-DD'
        ),

    # #### any covid test in comparison period
    # anytest_date=patients.with_test_result_in_sgss(
    #     pathogen="SARS-CoV-2",
    #     test_result="any",
    #     between=["start_k_date + 1 days", "end_k_date"],
    #     restrict_to_earliest_specimen_date=False,
    #     find_first_match_in_period=True,
    #     returning="date",
    #     date_format = "YYYY-MM-DD",
	#     ),


    #### at-risk group variables: updating
    # Asthma Admission codes
    astadm=patients.with_these_clinical_events(
      astadm_primis,
      returning="binary_flag",
      between=["start_k_date - 2 years", "start_k_date"],
    ),
    # Asthma systemic steroid prescription code in month 1
    astrxm1=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between=["start_k_date - 30 days", "start_k_date"],
    ),
    # Asthma systemic steroid prescription code in month 2
    astrxm2=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between=["start_k_date - 60 days", "start_k_date - 31 days"],
    ),
    # Asthma systemic steroid prescription code in month 3
    astrxm3=patients.with_these_medications(
      astrx_primis,
      returning="binary_flag",
      between=["start_k_date - 90 days", "start_k_date - 61 days"],
    ),

    # Chronic kidney disease diagnostic codes
    ckd_group=patients.satisfying(
        """
            ckd OR
            (ckd15_date AND 
            (ckd35_date >= ckd15_date) OR (ckd35_date AND NOT ckd15_date))
        """,
        # Chronic kidney disease codes - all stages
        ckd15_date=patients.with_these_clinical_events(
            ckd15_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="start_k_date",
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease codes-stages 3 - 5
        ckd35_date=patients.with_these_clinical_events(
            ckd35_primis,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="start_k_date",
            date_format="YYYY-MM-DD",
        ),
        # Chronic kidney disease diagnostic codes
        ckd=patients.with_these_clinical_events(
            ckd_primis,
            returning="binary_flag",
            on_or_before="start_k_date",
        ),
        return_expectations={"incidence": 0.01},
    ),

    # Immunosuppression medication codes
    immrx=patients.with_these_medications(
        immrx_primis,
        returning="binary_flag",
        between=["start_k_date - 180 days", "start_k_date"],
    ),

    # most recent numeric BMI
    bmi=patients.most_recent_bmi(
        on_or_before="start_k_date",
        minimum_age_at_measurement=16,
        include_measurement_date=True,
        date_format="YYYY-MM-DD",
        return_expectations={ "float": {"distribution": "normal", "mean": 28, "stddev": 8} },
        ),
    # most recent BMI stage
    bmi_stage=patients.with_these_clinical_events(
        bmi_stage_primis,
        on_or_before="start_k_date",
        returning="category",
        include_date_of_match=True,
        date_format="YYYY-MM-DD",
        find_last_match_in_period=True,
        return_expectations={
            "category": {"ratios": {"Obese class I": 0.4, "Obese class II": 0.4, "Obese class III": 0.2,}},
            "incidence": 0.75,
        },
        ),
    
    # pregnancy
    preg_group=patients.satisfying(
        """
        (preg_36wks_date) AND
        (pregdel_pre_date <= preg_36wks_date OR NOT pregdel_pre_date)
        """,
        # date of last pregnancy code in 36 weeks before ref_cev
        preg_36wks_date=patients.with_these_clinical_events(
            preg_primis,
            returning="date",
            find_last_match_in_period=True,
            between=["start_k_date - 252 days", "start_k_date"],
            date_format="YYYY-MM-DD",
        ),
        # date of last delivery code recorded in 36 weeks before elig_date
        pregdel_pre_date=patients.with_these_clinical_events(
            pregdel_primis,
            returning="date",
            find_last_match_in_period=True,
            between=["start_k_date - 252 days", "start_k_date"],
            date_format="YYYY-MM-DD",
        ),
    ),

    cev_group=patients.satisfying(
        "severely_clinically_vulnerable AND NOT less_vulnerable",

        # SHIELDED GROUP - first flag all patients with "high risk" codes
        severely_clinically_vulnerable=patients.with_these_clinical_events(
            shield_primis,
            returning="binary_flag",
            on_or_before="start_k_date",
            find_last_match_in_period=True,
        ),

        # find date at which the high risk code was added
        severely_clinically_vulnerable_date=patients.date_of(
            "severely_clinically_vulnerable",
            date_format="YYYY-MM-DD",
        ),

        # NOT SHIELDED GROUP (medium and low risk) - only flag if later than 'shielded'
        less_vulnerable=patients.with_these_clinical_events(
            nonshield_primis,
            between=["severely_clinically_vulnerable_date + 1 day", "start_k_date"],
        ),
        return_expectations={"incidence": 0.01},
    ),

#     housebound = patients.satisfying(
#     """
#     housebound_date
#     AND NOT no_longer_housebound
#     AND NOT moved_into_care_home
#     """,
        
#     housebound_date=patients.with_these_clinical_events( 
#       housebound, 
#       on_or_before="start_k_date",
#       find_last_match_in_period = True,
#       returning="date",
#       date_format="YYYY-MM-DD",
#     ),   
#     no_longer_housebound=patients.with_these_clinical_events( 
#       no_longer_housebound, 
#       on_or_after="housebound_date",
#     ),
#     moved_into_care_home=patients.with_these_clinical_events(
#       longres_primis,
#       on_or_after="housebound_date",
#     ),
#   ),

)
