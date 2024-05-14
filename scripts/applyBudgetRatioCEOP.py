"""Apply budget ratio to CEOP config"""
import argparse
import os
import numpy as np


def apply_budget_ratio_to_config(config_dir: str = "../configs/",
                                 config_fname: str = "ceop_config.json",
                                 budget_to_update: float = 0) -> None:
    """
    :param config_dir:
    :param config_fname:
    :param budget_ratio:
    :return:
    """
    with open(config_dir + config_fname, "r") as rf:
        lines = rf.read().splitlines()

    with open(config_dir + config_fname, "r+") as wf:
        for line in lines:
            if "budget" in line:
                wf.write("\t\t\"budget\": " + str(round(budget_to_update, 2)) + "\n")
            else:
                wf.write(line + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Apply budget ratio to CEOP config')
    parser.add_argument('--config_dir', type=str, default="../configs/")
    parser.add_argument('--config_fname', type=str, default="general.tddp.json")
    parser.add_argument('--instance', type=str, default="bubbles1")
    parser.add_argument('--budget_ratio', type=str, default="br100")
    args = parser.parse_args()

    fname_list = ["bonus1000", "bonus1000rdmRad", "bubbles1", "bubbles2", "bubbles3", "bubbles4", "bubbles5",
                  "bubbles6", "bubbles7", "bubbles8", "bubbles9", "chaoSingleDep", "concentricCircles1",
                  "concentricCircles2", "concentricCircles3", "concentricCircles4", "concentricCircles5", "d493rdmRad",
                  "dsj1000rdmRad", "kroD100rdmRad", "lin318rdmRad", "pcb442rdmRad", "rat195rdmRad", "rd400rdmRad",
                  "rotatingDiamonds1", "rotatingDiamonds2", "rotatingDiamonds3", "rotatingDiamonds4",
                  "rotatingDiamonds5", "team1_100", "team1_100rdmRad", "team2_200", "team2_200rdmRad", "team3_300",
                  "team3_300rdmRad", "team4_400", "team4_400rdmRad", "team5_499", "team5_499rdmRad", "team6_500",
                  "team6_500rdmRad", "d493or002", "dsj1000or002", "kroD100or002", "lin318or002", "pcb442or002",
                    "rat195or002", "rd400or002"]
    max_budget_list = [387.13, 932.07, 349.13, 428.28, 529.96, 805.46, 1038.16,
                       1229.66, 1607.31, 1946.72, 2259.22, 1039.61, 53.16,
                       153.13, 270.04, 452.64, 632.99, 134.28, 624.75, 141.83,
                       2047.42, 220.00, 68.22, 1246.69, 32.29, 140.48, 380.88,
                       770.66, 1510.75, 307.34, 388.54, 246.68, 613.66, 464.20,
                       378.09, 685.52, 1000.31, 700.50, 446.19, 225.22, 620.98,
                       202.23, 938.71, 159.04, 2838.54, 322.54, 158.32, 1032.04]

    if args.budget_ratio == "br120":
        bgt_rate = 1.2
    elif args.budget_ratio == "br90":
        bgt_rate = 0.9
    elif args.budget_ratio == "br60":
        bgt_rate = 0.6
    else:
        raise ValueError("Invalid budget ratio: " + args.budget_ratio)

    bgt2update = max_budget_list[fname_list.index(args.instance)] / 60 * bgt_rate
    apply_budget_ratio_to_config(args.config_dir, args.config_fname, bgt2update)
