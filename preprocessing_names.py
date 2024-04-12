# File used to preprocess names of VaMiGyn to include sample aliases instead of sequencing meta data
# Positivte controls should be uploaded but not negative controls

import pandas as pd


samples_s3_raw = pd.read_csv(
    "/ceph/projects/031_VaMiGyn/ENA_upload/host_removed/samples_all/all_samples.txt",
    header=None,
)

dysplasia_samples = pd.read_csv(
    "/ceph/projects/031_VaMiGyn/ENA_upload/List_dysplasia_samples.csv",
    header=0,
    sep=";",
)

control_samples = pd.read_csv(
    "/ceph/projects/031_VaMiGyn/ENA_upload/List_controls_samples.csv", header=0, sep=";"
)

print(samples_s3_raw)
print(dysplasia_samples)
