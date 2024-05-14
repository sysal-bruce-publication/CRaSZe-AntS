#!/usr/bin/env bash

# Declare problem name
PROBLEM="tddp"

# Declare instance name
CASE="bubbles1"

# Budget ratio option ["br120", "br90", "br60", "br30""] (note that "br30" is not available for TDDP)
BR="br60"

# Enter the directory of workspace
cd <path-to-workspace>

# Set budget for configuration json file
python3 <path-to-script>/applyBudgetRatioSOP.py --config_dir "configs/" --instance "$CASE.$PROBLEM" --budget_ratio "$BR"

# Argument for binary file: {json-file-name} {instance-name} {print configuration?} {write Steiner Zone vertex?}
workspace/crazyantsTDDP "example_config.$PROBLEM.json" "$CASE" 1 1