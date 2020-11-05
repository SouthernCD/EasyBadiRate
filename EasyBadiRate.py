import copy
import subprocess
import os
import re
import shutil
import logging
from Bio import Phylo
from io import StringIO

def logging_init(program_name, log_file=None):
    # create logger with 'program_name'
    logger = logging.getLogger(program_name)
    logger.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if not log_file is None:
        # create file handler which logs even debug messages
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)

    return logger

def rmdir(dir_name):
    if os.path.exists(dir_name):
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)
        else:
            os.remove(dir_name)

def mkdir(dir_name, keep=False):
    if keep is False:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        os.makedirs(dir_name)
    else:
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    return dir_name

def cmd_run(cmd_string, cwd=None, retry_max=5, silence=False, log_file=None):
    module_logger = logging_init("cmd_run", log_file)
    module_logger.info("Calling a bash cmd with retry_max %d: %s" %
                       (retry_max, cmd_string))
    if not silence:
        print("Running " + str(retry_max) + " " + cmd_string)
    p = subprocess.Popen(cmd_string, shell=True,
                         stderr=subprocess.PIPE, stdout=subprocess.PIPE, cwd=cwd)
    output, error = p.communicate()
    if not silence:
        print(error.decode())
    returncode = p.poll()
    module_logger.info("Finished bash cmd with returncode as %d" % returncode)
    if returncode == 1:
        if retry_max > 1:
            retry_max = retry_max - 1
            cmd_run(cmd_string, cwd=cwd, retry_max=retry_max)
    del module_logger.handlers[:]

    output = output.decode()
    error = error.decode()

    return (not returncode, output, error)

def get_sons(parent_clade):
    """
    get most close son nodes for clade in a tree
    """
    return parent_clade.clades

def run_badirate_one(ctl_file, badirate_path):
    cmd_run("perl %s %s" % (badirate_path, ctl_file), silence=True)


def run_badirate(ctl_list, output_list, badirate_path):
    # list_run=[ctl_file1,ctl_file2,..]
    for file_tmp in output_list:
        if os.path.exists(file_tmp):
            os.remove(file_tmp)

    for ctl_file in ctl_list:
        run_badirate_one(ctl_file, badirate_path)

    output_file_dict = {}

    for output_file in output_list:
        output_file_dict[output_file] = 0
        if os.path.exists(output_file):
            output_file_dict[output_file] = os.path.getsize(output_file)

    return output_file_dict


def get_K(bmodel_string):
    return 2*(len(re.findall(r'_', bmodel_string))+1)


def two_better_than_one_significance(k1, Likelihood1, k2, Likelihood2):
    AIC1 = 2*float(k1)-2*float(Likelihood1)
    AIC2 = 2*float(k2)-2*float(Likelihood2)

    if max([AIC1, AIC2])-min([AIC1, AIC2]) > 2:
        if AIC1 > AIC2:
            return 1, AIC1, AIC2
        elif AIC1 > AIC2:
            return -1, AIC1, AIC2
    else:
        return 0, AIC1, AIC2


def get_branch_name(badirate_tree, size_file, badirate_path):
    """
    perl ~/Program/badirate/badirate-master/BadiRate.pl -sizefile 1 -print_ids -treefile species_tree.tre > species_tree.badirate.tre
    """

    cmd_string = "perl %s -sizefile %s -print_ids -treefile %s" % (badirate_path, size_file, badirate_tree)
    flag, output, error = cmd_run(cmd_string,silence=True)

    labeled_tree_string = output.split('\n')[1]
    labeled_tree = StringIO(labeled_tree_string)

    tree = Phylo.read(labeled_tree, 'newick')

    all_branch = []

    for clade in tree.find_clades(order='level'):
        if clade.is_terminal():
            continue
        else:
            for son in get_sons(clade):
                clade_name = clade.confidence
                if son.is_terminal():
                    son_name = son.name.split("_")[-1]
                else:
                    son_name = son.confidence

                all_branch.append("%s->%s" % (clade_name, son_name))

    return all_branch, labeled_tree_string


def make_bmodel_string(all_branch, known="", new_ratio=0, GR=0, FR=0):
    """
    all_branch = ["17->15","15->11","11->9","9->3","3->1","3->2","9->8","8->6","6->4","6->5","8->7","11->10","15->14","14->12","14->13","17->16"]
    """
    all_string = {}

    if GR == 1:
        output = ""
        for i in all_branch:
            output = output+i+":"
        output = output.rstrip(":")
        return output

    if FR == 1:
        output = ""
        for i in all_branch:
            output = output+i+"_"
        output = output.rstrip("_")
        return output

    if known == "":
        head_string = known
        remain_branch = copy.deepcopy(all_branch)
    else:
        if new_ratio == 1:
            head_string = known+"_"
        elif new_ratio == 0:
            head_string = known+":"
        known_branchs = []
        ratios = known.split("_")
        for ratio in ratios:
            branchs = ratio.split(":")
            for branch in branchs:
                known_branchs.append(branch)

        remain_branch = copy.deepcopy(all_branch)
        for branch in known_branchs:
            remain_branch.remove(branch)

    output = {}
    for branch in remain_branch:
        output_string = head_string+branch+"_"
        for branch_other in remain_branch:
            if branch == branch_other:
                continue
            output_string = output_string + branch_other + ":"
        output_string = output_string.rstrip(":")
        output[branch] = output_string
    return head_string, output, new_ratio


def make_control_file(bmodel_string, size_file, treefile, work_dir, num):

    clt_0_file = "%s/%s.0.ctl" % (work_dir, str(num))
    out_0_file = "%s/%s.0.out" % (work_dir, str(num))

    with open(clt_0_file, 'w') as f:
        f.write(
            """
root_dist = 1
sizefile = %s
treefile = %s
n_max_int = 10
priorfile = 0
outlier = 0
seed = 587347092
unobs = 1
rmodel = GD
ep = ML
help = 0
out = %s
anc = 1
version = 0
print_ids = 0
bmodel = %s
start_val = 0
family = 0
""" % (size_file, treefile, out_0_file, bmodel_string))

    clt_1_file = "%s/%s.1.ctl" % (work_dir, str(num))
    out_1_file = "%s/%s.1.out" % (work_dir, str(num))

    with open(clt_1_file, 'w') as f:
        f.write(
            """
root_dist = 1
sizefile = %s
treefile = %s
n_max_int = 10
priorfile = 0
outlier = 0
seed = 1
unobs = 1
rmodel = GD
ep = ML
help = 0
out = %s
anc = 1
version = 0
print_ids = 0
bmodel = %s
start_val = 1
family = 0
""" % (size_file, treefile, out_1_file, bmodel_string))

    return [[clt_0_file, out_0_file], [clt_1_file, out_1_file]]


def badirate_output_parse(file_name):
    F1 = open(file_name)
    # print file_name
    all_text = F1.read()
    info = all_text.split('--------------------\n')
    # print info
    while '' in info:
        info.remove('')
    INTERNAL_ID_TREE = info[0]
    INPUT = info[1]
    OUTPUT = info[2]
    info_output = OUTPUT.split('\n')
    for i in info_output:
        match = re.match(r'\t\t#Likelihood: (\S+)', i)
        if match:
            Likelihood = match.group(1)
    # print Likelihood
    info_output = OUTPUT.split('##')
    flag = 0
    mini_dict = {}
    for i in info_output:
        match = re.match(r'Minimum number of gains and losses per branch', i)
        if match:
            min_num = i.split('\n')
            for line in min_num:
                match2 = re.match(r'\t\t(\d+\->\d+)\t(\d+)\t(\d+)', line)
                if match2:
                    mini_dict[match2.group(1)] = [int(
                        match2.group(2)), int(match2.group(3))]

    return INPUT, mini_dict, float(Likelihood)


def get_best_start_value(input_list):
    [like1, like2] = input_list
    if like1 == '-inf':
        like1 = -99999999999999999999
    else:
        like1 = float(like1)

    if like2 == '-inf':
        like2 = -99999999999999999999
    else:
        like2 = float(like2)

    temp_list = [like1, like2]
    best = max(temp_list)
    if best == -99999999999999999999:
        return "-inf"
    else:
        return temp_list.index(best)


def main_pipeline(tag, size_file, species_tree, work_dir, badirate_path, label_tree_path, keep_tmp_dir):
    log_file = work_dir + "/log"

    mkdir(work_dir, False)

    logger = logging_init("badirate_exp_con", log_file)
    logger.info("Build work dir")

    all_branch, label_tree = get_branch_name(species_tree, size_file, badirate_path)

    if not label_tree_path is None:
        with open(label_tree_path, 'w') as f:
            f.write(label_tree)

    # cal free branch model
    logger.info("cal free branch model")

    FR_string = make_bmodel_string(all_branch, FR=1)
    out_list = make_control_file(
        FR_string, size_file, species_tree, work_dir, "FR")

    a = [out_list[0][0], out_list[1][0]]
    b = [out_list[0][1], out_list[1][1]]

    output_file_dict = run_badirate(a, b, badirate_path)
    logger.info("cal free branch model, done!")

    output_file_list = list(output_file_dict)

    like_FR = [badirate_output_parse(file)[2] for file in output_file_list]

    best_like_FR_index = get_best_start_value(like_FR)

    if best_like_FR_index == "-inf":
        logger.info("%s can't get information from free model, failed!" % tag)
        return None
    else:
        best_like_FR = like_FR[best_like_FR_index]

    FR_output_file = output_file_list[best_like_FR_index]

    INPUT, mini_dict, Likelihood = badirate_output_parse(FR_output_file)
    back_branch = []
    test_branch = []
    for i in mini_dict:
        gains, losses = mini_dict[i]
        if gains == losses and gains == 0:
            back_branch.append(i)
        else:
            test_branch.append(i)

    logger.info("free branch model get likelihood: %.5f, and %d branch have change" % (
        Likelihood, len(all_branch) - len(back_branch)))
    # not branch changed
    if len(back_branch) == len(all_branch):
        logger.info("No exp & con")
        return [tag,[],[],Likelihood]

    like_FR = sorted(like_FR, reverse=True)

    logger.info("cal Null hypothesis branch model")

    eFR_model_string = ""
    for i in back_branch:
        eFR_model_string = eFR_model_string+i+":"
    eFR_model_string = eFR_model_string.rstrip(":")
    back_branch_string = eFR_model_string
    eFR_model_string = eFR_model_string+"_"

    for i in all_branch:
        if i in back_branch:
            continue
        eFR_model_string = eFR_model_string+i+"_"
    eFR_model_string = eFR_model_string.rstrip("_")

    logger.info("get eFR model string: %s" % eFR_model_string)
    out_list = make_control_file(
        eFR_model_string, size_file, species_tree, work_dir, "eFR")

    a = [out_list[0][0], out_list[1][0]]
    b = [out_list[0][1], out_list[1][1]]

    output_file_dict = run_badirate(a, b, badirate_path)
    logger.info("cal Null hypothesis branch model, done!")

    output_file_list = list(output_file_dict)

    like_eFR = [badirate_output_parse(file)[2] for file in output_file_list]

    best_like_eFR_index = get_best_start_value(like_eFR)
    if best_like_eFR_index == "-inf":
        logger.info(
            "%s can't get information from null hypothesis free model, failed!" % tag)
        return None
    else:
        best_like_eFR = like_eFR[best_like_eFR_index]

    eFR_output_file = output_file_list[get_best_start_value(like_eFR)]

    INPUT, mini_dict, Likelihood = badirate_output_parse(eFR_output_file)

    output_dict = {}
    output_dict["family_id"] = tag
    output_dict["FR"] = FR_output_file
    output_dict["eFR"] = {}
    output_dict["tsv_file"] = size_file
    output_dict["eFR"]['Likelihood'] = Likelihood
    output_dict["eFR"]['K'] = get_K(eFR_model_string)
    output_dict["eFR"]['output_file'] = eFR_output_file

    logger.info("null hypothesis branch model get likelihood: %.5f, K is %d" % (
        Likelihood, output_dict["eFR"]['K']))

    # begin test

    good_branch = []
    for branch_now in test_branch:

        logger.info("test %s %s now:" % (branch_now, mini_dict[branch_now]))

        logger.info("cal test branch model")
        num = re.sub("->", '_', branch_now)
        output_dict[num] = {}
        output_dict[num]['branch'] = branch_now
        test_branch_string = back_branch_string+":"+branch_now+"_"
        for j in test_branch:
            if not j == branch_now:
                test_branch_string = test_branch_string+j+"_"
        test_branch_string = test_branch_string.rstrip("_")

        logger.info("get test branch model string: %s" % test_branch_string)

        out_list = make_control_file(
            test_branch_string, size_file, species_tree, work_dir, num)

        a = [out_list[0][0], out_list[1][0]]
        b = [out_list[0][1], out_list[1][1]]

        output_file_dict = run_badirate(a, b, badirate_path)
        logger.info("cal test branch model, done!")

        like_now = [(badirate_output_parse(file[1])[2], file)
                    for file in out_list]
        like_now = sorted(like_now, reverse=True)
        best_eFR = like_now[0]

        if best_eFR[0] == '-inf':
            output_dict[num]['out_file'] = best_eFR[1][1]
            output_dict[num]['Likelihood'] = '-99999999999999999999'
            output_dict[num]['bmodel_string'] = test_branch_string
            output_dict[num]['K'] = get_K(test_branch_string)
        else:
            output_dict[num]['out_file'] = best_eFR[1][1]
            output_dict[num]['Likelihood'] = best_eFR[0]
            output_dict[num]['bmodel_string'] = test_branch_string
            output_dict[num]['K'] = get_K(test_branch_string)

        logger.info("comp %s (%.5f, %d) and eFR (%.5f, %d)" % (
            branch_now, output_dict[num]['Likelihood'], output_dict[num]['K'], output_dict["eFR"]['Likelihood'], output_dict["eFR"]['K']))

        flag, AIC1, AIC2 = two_better_than_one_significance(
            output_dict[num]['K'], output_dict[num]['Likelihood'], output_dict["eFR"]['K'], output_dict["eFR"]['Likelihood'])

        logger.info("AIC: %.5f vs %.5f flag: %d" % (AIC1, AIC2, flag))

        if flag == 1:
            output_dict[num]['significative'] = 1
            good_branch.append(branch_now)
        else:
            output_dict[num]['significative'] = 0

    up_down={}
    up_down['up']=[]
    up_down['down']=[]
    for i in good_branch:
        if mini_dict[i][0] > 0:
            up_down['up'].append(i)
        elif mini_dict[i][1] >0:
            up_down['down'].append(i)

    printer=""
    for i in up_down['up']:
        printer=printer+i+","
    printer = printer.rstrip(",")
    printer = printer+"\t"
    for i in up_down['down']:
        printer=printer+i+","
    printer = printer.rstrip(",")
    printer = printer + "\t" + str(output_dict["eFR"]['Likelihood'])
    printer = tag +"\t"+ printer


    print("Tag\tGain\tLoss\tLikelihood")
    print(tag,up_down['up'],up_down['down'],output_dict["eFR"]['Likelihood'])

    if keep_tmp_dir:
        print("temp_dir is : %s" % work_dir)
    else:
        rmdir(work_dir)


if __name__ == "__main__":
    import argparse
    import uuid


    parser = argparse.ArgumentParser(
        prog='EasyBadiRate', description='EasyBadiRate\n'
    )

    parser.add_argument('tag', type=str, help='a tag for this family')
    parser.add_argument('tree_file', type=str, help='a species tree file in newick')
    parser.add_argument('size_tsv_file', type=str, help='a gene family size file in tsv file')
    parser.add_argument('-l', '--label_tree', type=str, help='output a labeled tree from BadiRate', default=None)
    parser.add_argument('-k', '--keep_tmp_dir', help='keep temp running dir', action='store_true')

    args = parser.parse_args()

    script_dir_path = os.path.split(os.path.realpath(__file__))[0]
    badirate_path = script_dir_path + "/badirate/BadiRate.pl"
    main_pipeline(args.tag, args.size_tsv_file, args.tree_file, "/tmp/" + uuid.uuid4().hex, badirate_path, args.label_tree, args.keep_tmp_dir)
