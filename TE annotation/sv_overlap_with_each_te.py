import os
import sys
import numpy as np
import subprocess


intersect_output = sys.argv[1]
final_output = sys.argv[2]
ratio_of_SV =float(sys.argv[3])

out = open(final_output, "w")

temp_file = os.path.dirname(intersect_output) + "/temp.bed"

def run_bedtools_merge(bed_file):
    cmd ='sortBed -i ' + bed_file + ' | bedtools merge -i - -c 4 -o collapse > ' + bed_file + '.merged'
    os.system(cmd)
    output = [line.rstrip() for line in open(bed_file + '.merged').readlines()]
    s_len = 0
    all_tes_len = []
    all_tes_type = []
    for o in output:
        items = o.split("\t")
        s_len = s_len + int(items[2]) -int(items[1])
        all_tes_len += [int(i.split(":")[1]) for i in items[3].split(",")]
        all_tes_type += [i.split(":")[0] for i in items[3].split(",")]
    inx = np.argmax(np.array(all_tes_len))
    best_te = all_tes_type[inx]
    return best_te, s_len


def merge_ranges(bed_list):
    sort_bed_list = list(sorted(bed_list))
    # result_list = []
    result_len= 0
    low_est = sort_bed_list[0][0]
    high_est = sort_bed_list[0][1]
    for i in range(1, len(sort_bed_list)):
        low, high = sort_bed_list[i]
        if high_est >= low:
            if high_est < high:
                high_est = high
        else:
            # result_list.append([low_est, high_est])
            result_len += high_est-low_est
            low_est, high_est = sort_bed_list[i]
    # result_list.append([low_est, high_est])
    result_len += high_est-low_est
    return result_len

# chr3	31320567	31334095	chr3	31320557	31321002	Others	435
# chr3	31320567	31334095	chr3	31321032	31321126	Others	94
# chr3	31320567	31334095	chr3	31321065	31321820	Gypsy	755
# chr3	31320567	31334095	chr3	31321606	31323235	Others	1629
# chr3	31320567	31334095	chr3	31330935	31331663	Gypsy	728
# chr3	31320567	31334095	chr3	31331664	31332500	Gypsy	836
# chr3	31320567	31334095	chr3	31331664	31333065	Others	1401

data = [line.rstrip().split("\t") for line in open(intersect_output, "r").readlines()]

sv_te ={}
for d in data:
    temp_key = d[0] + "\t" + d[1] + "\t" + d[2]
    if temp_key in sv_te:
        if d[6] in sv_te[temp_key]:
            sv_te[temp_key][d[6]] += [[max(int(d[1]), int(d[4])), min(int(d[2]), int(d[5]))]]
        else:  
            sv_te[temp_key][d[6]] = [[max(int(d[1]), int(d[4])), min(int(d[2]), int(d[5]))]]
    else:
        sv_te[temp_key] = {}
        sv_te[temp_key][d[6]] = [[max(int(d[1]), int(d[4])), min(int(d[2]), int(d[5]))]]

for current_key in sv_te:
    sv_len = (int(current_key.split("\t")[2]) - int(current_key.split("\t")[1]))*ratio_of_SV
    out.write(current_key + "\t")
    for each_te in sv_te[current_key]:
        current_len = merge_ranges(sv_te[current_key][each_te])
        if current_len > sv_len:
            out.write(each_te + ":" + str(current_len) + "|")
    out.write("\n")
out.close()