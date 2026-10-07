import re


include: "constants.smk"
include: "utils.smk"


configfile: "config/config.yaml"


pepfile: "config/pep/config.yaml"


if config["dir_run"] and config["dir_run"] is not None:

    workdir: config["dir_run"]


DIR_DATA = config["dir_data"]
DF_SAMPLE = pep.sample_table.set_index("sample_name", drop=False)
if "library_layout" not in DF_SAMPLE.columns:
    DF_SAMPLE["library_layout"] = "paired-end"

SAMPLES = DF_SAMPLE.index
SAMPLES_PE = SAMPLES[DF_SAMPLE["library_layout"] == "paired-end"]
SAMPLES_SE = SAMPLES[DF_SAMPLE["library_layout"] == "single-end"]

TO_CLEAN_FQ = config["clean_fq"]
TO_RUN_FASTQC = config["run_fastqc"]
TO_RUN_MULTIQC = config["run_multiqc"]


def _wildcard_constraint(samples):
    samples = [str(s) for s in samples]
    return r"|".join(map(re.escape, samples)) if samples else r"(?!)"


wildcard_constraints:
    sample=_wildcard_constraint(SAMPLES),
    sample_pe=_wildcard_constraint(SAMPLES_PE),
    sample_se=_wildcard_constraint(SAMPLES_SE),
    type_fragment=r"(nfr|1n|2n|3n)",
    caller=r"(macs3|genrich)",
