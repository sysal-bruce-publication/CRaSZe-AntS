import numpy as np
import argparse


def change_radius_prize(data_dir: str, raw_fname: str, r: float=-1) -> None:
    """
    Change the radii to fixed value and restrict the prize to [0.8, 1.0].
    If unit prize, then set all to 0.9.
    :param data_dir:
    :param raw_fname:
    :param r:
    :return:
    """
    str_xs, str_ys, not_used1, raw_rs, raw_ps = [], [], [], [], []
    with open(data_dir + raw_fname, 'r') as rf:
        lines = rf.read().splitlines()
        for line in lines:
            line = line.split(" ")
            str_xs.append(line[0])
            str_ys.append(line[1])
            not_used1.append(line[2])
            raw_rs.append(float(line[3]))
            raw_ps.append(float(line[4]))
    # Change radius first
    raw_rs = np.array(raw_rs, dtype=np.float64)
    if r != -1:
        raw_rs.fill(r)
    raw_rs[0] = raw_rs[1] = 0

    # Change prize then
    raw_ps[0] = raw_ps[1] = 0
    # Normalize the prize list
    raw_ps = np.array(raw_ps, dtype=np.float64)
    p_max, p_min = np.max(raw_ps[2:]), np.min(raw_ps[2:])
    if p_max == p_min:
        for i in range(2, len(raw_ps)):
            raw_ps[i] = 1
    else:
        for i in range(2, len(raw_ps)):
            raw_ps[i] = (raw_ps[i] - p_min) / (p_max - p_min)

    # Change weights then
    weights = np.zeros_like(raw_ps, dtype=np.float64)

    w_max, w_min = np.max(weights[2:]), np.min(weights[2:])
    if w_max == w_min:
        for i in range(2, len(weights)):
            weights[i] = 0.9
    else:
        W_MAX, W_MIN = 1.0, 0.8
        for i in range(2, len(weights)):
            weights[i] = np.random.uniform(W_MIN, W_MAX)

    with open(data_dir + raw_fname, "w") as wf:
        for i in range(len(str_xs)):
            wf.write(str_xs[i] + " " + str_ys[i] + " " + not_used1[i]
                     + " " + str(raw_rs[i]) + " " + str(round(raw_ps[i], 5))
                     + " " + str(round(weights[i], 5)) + "\n")


def apply_ovlp_ratio_to_instances(data_dir: str, raw_fname: str, ovlp_ratio=0.1) -> None:
    """

    :param data_dir:
    :param raw_fname:
    :param ovlp_ratio:
    :return:
    """
    str_xs, str_ys, raw_xs, raw_ys, not_used1, raw_rs, not_used2 = [], [], [], [], [], [], []
    with open(data_dir + raw_fname, 'r') as rf:
        lines = rf.read().splitlines()
        for line in lines:
            line = line.split(" ")
            str_xs.append(line[0])
            str_ys.append(line[1])
            raw_xs.append(float(line[0]))
            raw_ys.append(float(line[1]))
            not_used1.append(line[2])
            raw_rs.append(line[3])
            not_used2.append(line[4])
    raw_xs = np.array(raw_xs, dtype=np.float64)
    raw_ys = np.array(raw_ys, dtype=np.float64)
    raw_rs = np.array(raw_rs, dtype=np.float64)
    x_max, x_min = np.max(raw_xs), np.min(raw_xs)
    y_max, y_min = np.max(raw_ys), np.min(raw_ys)
    l_contain = max(x_max - x_min, y_max - y_min)
    r = ovlp_ratio * l_contain  # / (1 - 2 * ovlp_ratio)
    raw_rs.fill(r)
    raw_rs[0] = 0.

    instance_name = raw_fname.split(".")[0]
    problem_name = raw_fname.split(".")[1]
    or_str = instance_name + "or" + str(ovlp_ratio).replace('.', '') + "." + problem_name
    with open(data_dir + or_str, "w") as wf:
        for i in range(len(str_xs)):
            wf.write(str_xs[i] + " " + str_ys[i] + " " + not_used1[i]
                     + " " + str(raw_rs[i]) + " " + not_used2[i] + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Solution Plot')
    parser.add_argument('--dir', type=str, default="../data/ceop/instances/")
    parser.add_argument('-f', '--filename', type=str, default="bubbles4.ceop")
    parser.add_argument('--overlap', type=str, default="or01")
    parser.add_argument('-r', '--radius', type=float, default=5)
    args = parser.parse_args()

    if args.radius <= 0:
        if args.overlap == "or001":
            ovlp = 0.01
        elif args.overlap == "or002":
            ovlp = 0.02
        elif args.overlap == "or003":
            ovlp = 0.03
        elif args.overlap == "or004":
            ovlp = 0.04
        elif args.overlap == "or005":
            ovlp = 0.05
        elif args.overlap == "or01":
            ovlp = 0.1
        elif args.overlap == "or03":
            ovlp = 0.3
        else:
            raise Exception("Invalid overlap ratio")
        apply_ovlp_ratio_to_instances(data_dir=args.dir,
                                      raw_fname=args.filename, ovlp_ratio=ovlp)
    else:
        change_radius_prize(data_dir=args.dir, raw_fname=args.filename,
                            r=args.radius)

