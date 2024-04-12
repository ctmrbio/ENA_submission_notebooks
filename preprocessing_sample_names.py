# File used to preprocess names of VaMiGyn to include sample aliases instead of sequencing meta data
# Positivte controls should be uploaded but not negative controls

import pandas as pd
import numpy as np

# Load data
samples_s3_raw = pd.read_csv(
    "/ceph/projects/031_VaMiGyn/ENA_upload/host_removed/samples_all/all_samples.txt",
    header=0,
)

dysplasia_samples = pd.read_csv(
    "/ceph/projects/031_VaMiGyn/ENA_upload/List_dysplasia_samples.csv",
    header=0,
    sep=";",
)

# Add a tailing "__" for easier sub-string matching
dysplasia_samples["Flowcell_w_us"] = dysplasia_samples["Flowcell_L1_bc"] + "__"

control_samples = pd.read_csv(
    "/ceph/projects/031_VaMiGyn/ENA_upload/List_controls_samples.csv", header=0, sep=";"
)

# Add a tailing "__" for easier sub-string matching
control_samples["Flowcell_w_us"] = control_samples["Flowcell_L1_bc"] + "__"

# Seperate positive controls, negative controls and samples from each other
positive_controls = samples_s3_raw[samples_s3_raw["Samples_s3"].str.contains("pos")]
negative_controls = samples_s3_raw[samples_s3_raw["Samples_s3"].str.contains("neg")]
samples_raw = samples_s3_raw[
    ~samples_s3_raw["Samples_s3"].str.contains("|".join(["pos", "neg", "empty", "ext"]))
]

# Subset dysplasia samples from S3 raw output
dysplasia_processed = samples_raw[
    samples_raw["Samples_s3"].apply(
        lambda x: any(
            substring in x for substring in dysplasia_samples["Flowcell_w_us"]
        )
    )
]

dysplasia_processed["unique_samples_names"] = dysplasia_processed[
    "Samples_s3"
].str.replace(r".{7}$", "")

dysplasia_processed["Flowcell_w_us"] = dysplasia_processed[
    "unique_samples_names"
].str.replace(r"^.{25}", "")

dysplasia_intermediary = pd.DataFrame(
    {"Samples_s3": [None] * 177, "Flowcell_w_us": [None] * 177}
)

# Create an intermediary file where all unique enteries are saved, to be later merged with control samples
dysplasia_intermediary["Samples_s3"] = dysplasia_processed[
    "unique_samples_names"
].unique()
dysplasia_intermediary["Flowcell_w_us"] = dysplasia_processed["Flowcell_w_us"].unique()

dysplasia_merged = pd.merge(
    dysplasia_samples, dysplasia_intermediary, on="Flowcell_w_us", how="left"
)

# Complie a new name where the leading string in this name is the sampleID
dysplasia_merged["new_name"] = (
    dysplasia_merged["SampleID"] + "_" + dysplasia_merged["Flowcell_w_us"]
)

#######################################################################################
# Subset control samples from S3 raw output
control_processed = samples_raw[
    samples_raw["Samples_s3"].apply(
        lambda x: any(substring in x for substring in control_samples["Flowcell_w_us"])
    )
]

control_processed["unique_samples_names"] = control_processed["Samples_s3"].str.replace(
    r".{7}$", ""
)

control_processed["Flowcell_w_us"] = control_processed[
    "unique_samples_names"
].str.replace(r"^.{25}", "")

control_intermediary = pd.DataFrame(
    {"Samples_s3": [None] * 177, "Flowcell_w_us": [None] * 177}
)

# Create an intermediary file where all unique enteries are saved, to be later merged with control samples
control_intermediary["Samples_s3"] = control_processed["unique_samples_names"].unique()
control_intermediary["Flowcell_w_us"] = control_processed["Flowcell_w_us"].unique()

control_merged = pd.merge(
    control_samples, control_intermediary, on="Flowcell_w_us", how="left"
)

# Complie a new name where the leading string in this name is the sampleID
control_merged["new_name"] = (
    control_merged["SampleID"] + "_" + control_merged["Flowcell_w_us"]
)

#######################################################################################
dysplasia_new_names_doubled = pd.DataFrame(np.repeat(dysplasia_merged["new_name"], 2))
dysplasia_new_names_doubled.index = range(len(dysplasia_new_names_doubled))
dysplasia_new_names_doubled["Flowcell_w_us"] = dysplasia_new_names_doubled[
    "new_name"
].str.replace(r"(^.{10})", "")
dysplasia_new_names_doubled = dysplasia_new_names_doubled[["Flowcell_w_us", "new_name"]]

dysplasia_processed.index = range(len(dysplasia_processed))

dysplasia_processed_merged = pd.merge(
    dysplasia_processed, dysplasia_new_names_doubled, on="Flowcell_w_us", how="left"
)

dysplasia_processed_merged = dysplasia_processed_merged.drop_duplicates()
dysplasia_processed.index = range(len(dysplasia_processed))

# Add tailing sequencing indecies
for value in range(len(dysplasia_processed_merged)):
    if value % 2 == 0:
        dysplasia_processed_merged.iloc[
            value, dysplasia_processed_merged.columns.get_loc("new_name")
        ] = (
            dysplasia_processed_merged.iloc[
                value, dysplasia_processed_merged.columns.get_loc("new_name")
            ]
            + "1.fq.gz"
        )
    else:
        dysplasia_processed_merged.iloc[
            value, dysplasia_processed_merged.columns.get_loc("new_name")
        ] = (
            dysplasia_processed_merged.iloc[
                value, dysplasia_processed_merged.columns.get_loc("new_name")
            ]
            + "2.fq.gz"
        )


#######################################################################################
control_new_names_doubled = pd.DataFrame(np.repeat(control_merged["new_name"], 2))
control_new_names_doubled.index = range(len(control_new_names_doubled))

print(control_merged)
control_new_names_doubled["Flowcell_w_us"] = control_new_names_doubled[
    "new_name"
].str.replace(r"(^.{9})", "")
control_new_names_doubled = control_new_names_doubled[["Flowcell_w_us", "new_name"]]
print(control_new_names_doubled)


# Add tailing sequencing indecies
# for value in range(len(control_processed_merged)):
#    if value % 2 == 0:
#        control_processed_merged.iloc[
#            value, control_processed_merged.columns.get_loc("new_name")
#        ] = (
#            control_processed_merged.iloc[
#                value, control_processed_merged.columns.get_loc("new_name")
#            ]
#            + "1.fq.gz"
#        )
#    else:
#        control_processed_merged.iloc[
#            value, control_processed_merged.columns.get_loc("new_name")
#        ] = (
#            control_processed_merged.iloc[
#                value, control_processed_merged.columns.get_loc("new_name")
#            ]
#            + "2.fq.gz"
#        )
