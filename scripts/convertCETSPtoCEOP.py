"""Convert CETSP instances to CEOP instances"""
import os


def convert_cetsp_to_ceop(inst_dir: str, inst_name: str):
    """
    :param inst_dir:
    :param inst_name:
    :return:
    """
    xs, ys, zs, rs, ps = [], [], [], [], []
    with open(inst_dir + inst_name + ".cetsp", 'r') as rf:
        lines = rf.read().splitlines()
        for line in lines:
            line = line.split(" ")
            xs.append(line[0])
            ys.append(line[1])
            zs.append(line[2])
            rs.append(line[3])
            ps.append(line[4])

    start_node = [xs[0], ys[0], zs[0], rs[0], ps[0]]
    xs.insert(0, start_node[0])
    ys.insert(0, start_node[1])
    zs.insert(0, start_node[2])
    rs.insert(0, start_node[3])
    ps.insert(0, start_node[4])
    with open(inst_dir + inst_name + ".ceop", "w") as wf:
        for i in range(len(xs)):
            wf.write(xs[i] + " " + ys[i] + " " + zs[i] + " "
                     + rs[i] + " " + ps[i] + "\n")


if __name__ == "__main__":
    fname_list = ["bonus1000", "bonus1000rdmRad", "bubbles1", "bubbles2", "bubbles3", "bubbles4", "bubbles5",
                  "bubbles6", "bubbles7", "bubbles8", "bubbles9", "chaoSingleDep", "concentricCircles1",
                  "concentricCircles2", "concentricCircles3", "concentricCircles4", "concentricCircles5", "d493rdmRad",
                  "dsj1000rdmRad", "kroD100rdmRad", "lin318rdmRad", "pcb442rdmRad", "rat195rdmRad", "rd400rdmRad",
                  "rotatingDiamonds1", "rotatingDiamonds2", "rotatingDiamonds3", "rotatingDiamonds4",
                  "rotatingDiamonds5", "team1_100", "team1_100rdmRad", "team2_200", "team2_200rdmRad", "team3_300",
                  "team3_300rdmRad", "team4_400", "team4_400rdmRad", "team5_499", "team5_499rdmRad", "team6_500",
                  "team6_500rdmRad"]
    suffix = ".cetsp"
    for fname in fname_list:
        if os.path.isfile("./" + fname + suffix):
            convert_cetsp_to_ceop("../instances/CEOP/", fname)

    ovlp_fname_list = ["d493", "dsj1000", "kroD100", "lin318", "pcb442", "rat195", "rd400"]
    ovlp_list = ["or03", "or01", "or002"]
    for fname in ovlp_fname_list:
        for ovlp in ovlp_list:
            if os.path.isfile("./" + fname + ovlp + suffix):
                convert_cetsp_to_ceop("../instances/CEOP/", fname + ovlp)