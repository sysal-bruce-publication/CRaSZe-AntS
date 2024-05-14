"""Convert CETSP instances to CEOP instances"""
import numpy as np
import os
import argparse


def convert_ceop_to_tdd(inst_dir: str, inst_name: str, r: float = -1,
                        min_p: float = 0.1, max_p: float = 0.9,
                        min_weight: float = 0.8):
    """
    :param max_p:
    :type max_p:
    :param min_p:
    :type min_p:
    :param min_weight:
    :type min_weight:
    :param r:
    :type r:
    :param inst_dir:
    :param inst_name:
    :return:
    """
    str_xs, str_ys, not_used1, raw_rs, raw_ps = [], [], [], [], []
    with open(inst_dir + inst_name + ".ceop", 'r') as rf:
        lines = rf.read().splitlines()
        for line in lines:
            line = line.split(" ")
            str_xs.append(line[0])
            str_ys.append(line[1])
            not_used1.append(line[2])
            raw_rs.append(float(line[3]))
            raw_ps.append(float(line[4]))
    # Radius
    raw_rs = np.array(raw_rs, dtype=np.float64)
    if r > 0:
        raw_rs.fill(r)
    raw_rs[0] = raw_rs[1] = 0
    # Prize
    raw_ps = np.array(raw_ps, dtype=np.float64)
    raw_p_max, raw_p_min = np.max(raw_ps[2:]), np.min(raw_ps[2:])
    if raw_p_max == raw_p_min:
        for i in range(2, len(raw_ps)):
            raw_ps[i] = 0.9
    else:
        factor = (max_p - min_p) / (raw_p_max - raw_p_min)
        for i in range(2, len(raw_ps)):
            raw_ps[i] = (raw_ps[i] - raw_p_min) * factor + min_p
    # Weight
    weights = np.zeros_like(raw_ps, dtype=np.float64)
    for i in range(2, len(weights)):
        weights[i] = np.random.uniform(min_weight, 1)

    with open(inst_dir + inst_name + ".tddp", "w") as wf:
        for i in range(len(str_xs)):
            wf.write(str_xs[i] + " " + str_ys[i] + " " + not_used1[i] + " "
                     + str(round(raw_rs[i], 5)) + " "
                     + str(round(raw_ps[i], 5)) + " "
                     + str(round(weights[i], 5)) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert CEOP instances to TDD instances')
    parser.add_argument('--inst_dir', type=str, default="./")
    parser.add_argument('--inst_name', type=str, default="team1_100")
    parser.add_argument('-r', '--radius', type=float, default=-1)
    args = parser.parse_args()

    fname_list = ["bonus1000", "chaoSingleDep",

                  "d493or002", "d493or01", "d493or03",
                  "dsj1000or002", "dsj1000or01", "dsj1000or03",
                  "kroD100or002", "kroD100or01", "kroD100or03",
                  "lin318or002", "lin318or01", "lin318or03",
                  "pcb442or002", "pcb442or01", "pcb442or03",
                  "rat195or002", "rat195or01", "rat195or03",
                  "rd400or002", "rd400or01", "rd400or03",

                  "bubbles1", "bubbles2", "bubbles3", "bubbles4", "bubbles5",
                  "bubbles6", "bubbles7", "bubbles8", "bubbles9",

                  "concentricCircles1", "concentricCircles2", "concentricCircles3",
                  "concentricCircles4", "concentricCircles5",

                  "rotatingDiamonds1", "rotatingDiamonds2", "rotatingDiamonds3",
                  "rotatingDiamonds4", "rotatingDiamonds5",

                  "team1_100", "team2_200", "team3_300",
                  "team4_400", "team5_499", "team6_500"]

    if args.inst_name in fname_list:
        if os.path.isfile("./" + args.inst_name + ".ceop"):
            convert_ceop_to_tdd(args.inst_dir, args.inst_name, args.radius)
