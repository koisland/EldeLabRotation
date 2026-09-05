from os.path import join, dirname

SAMPLES = config["samples"]
SAMPLE_NAMES = [sm_info["name"] for sm_info in SAMPLES]


wildcard_constraints:
    sm="|".join(SAMPLE_NAMES),
