#!/usr/bin/env python3
# Chongwei 20200601
# Email: chongwei.bi@kaust.edu.sa

import os
import sys
import time
import subprocess
import glob
import shutil
import argparse
import re
import logging
import multiprocessing
import pandas as pd
from _version import version
from datetime import datetime
from timeit import default_timer as timer
import call_variant


def get_argparse():
    parser = argparse.ArgumentParser(description='Real-Time analysis of Nanopore data for Covid-19 sequencing.')
    parser.add_argument('-p', '--path', type=str, required=True, help='path/to/nanopore_result_folder')
    parser.add_argument('-s', '--save_path', type=str,
                        help='path/to/saved_folder Default: rtnano_result in -p PATH folder')
    parser.add_argument('-r', '--refer_seq', type=str,
                        help='path/to/reference_genome.fa, default is using SARS-CoV-2.fa in program folder')
    parser.add_argument('-t', '--thread', type=int, default='1', help='working thread [1]')
    parser.add_argument('-T', '--interval_time', type=int, default='1',
                        help='interval time for analysis in minutes [1]')
    parser.add_argument('-R', '--target_region', type=str, default='11050-11244,14307-14500,23123-23431,28086-28752',
                        help='primer targeted region (comma separated). Default: '
                             '11050-11244,14307-14500,23123-23431,28086-28752')
    # 11075-11221,14329-14477,23144-23411,28111-28731

    parser.add_argument('-g', '--guppy_barcoder', type=str,
                        help='Optional: path/to/guppy_barcoder, when offering this parameter, it will do additional '
                             'demultiplexing using guppy_barcoder --require_barcodes_both_ends')
    parser.add_argument('--run_time', type=int, default='48',
                        help='total run time in hours [48]')
    parser.add_argument('--resume', action='store_true',
                        help='resume the unexpectedly interrupted analysis. Please use the same [-p] [-s] as before')
    parser.add_argument('--put_back', action='store_true',
                        help='return fastq file to their original fastq_pass folder. Please use the same [-p] [-s] [-g]'
                             ' as you generated the result, together with --put_back')
    parser.add_argument('--call_variant', action='store_true',
                        help='call variants using samtools and filter by alleic frequency (>=0.5). Please use the same [-p] [-s] '
                             'as you generated the result. It uses fastq file in analyzed_achieve/accumulated_reads '
                             'folder. If you want to use it for your own data, '
                             'please put fastq file in this folder, one sample one fastq file')
    parser.add_argument('--alle_freq', type=float, default='0.5', help='filter SNVs by allelic frequency [0.5]')
    parser.add_argument('-v', '--version', action='version', version=version, help='show the current version')
    args = parser.parse_args()
    if args.refer_seq is None:
        args.refer_seq = sys.path[0] + '/SARS-CoV-2.fa'

    return args


def prepare_env(args):
    input_path = args.path

    if args.save_path is None:
        result_folder = input_path + '/rtnano_result'
    else:
        result_folder = args.save_path

    if os.path.isdir(result_folder):
        if args.resume is not True:
            sys.stderr.write("WARNING!  Result folder exist! <%s> \n"
                             "WARNING!  To avoid overwriting in existing data, Please use a new path and restart the program.\n"
                             "WARNING!  If you want to resume a unexpectedly interrupted analysis, "
                             "please use the same save_path [-s] as before, together with --resume \n" % result_folder)
            sys.exit(1)
        else:
            return result_folder

    else:
        os.mkdir(result_folder)
        return result_folder


def get_fastq_file(args, result_folder):
    ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
    fastq_path = args.path + '/fastq_pass'

    analyzing = result_folder + '/analyzing/' + ctime
    analyzed = result_folder + '/analyzed_achieve'

    accumulated_reads = result_folder + '/analyzed_achieve/accumulated_reads'

    if not os.path.isdir(analyzed):
        os.mkdir(analyzed)
        pooled_result = analyzed + '/pooled_result.txt'
        with open(pooled_result, 'w') as outfile:
            header = ["#barcode", "%_mapped_record", "%_mapped_base", "%_ontarget_base", "read_number", "total_base",
                      "mapped_record", "unmapped_record", "mapped_base", "ontarget_base"]

            target_region = args.target_region.split(',')
            sum_header = header + target_region

            outfile.writelines('\t'.join(sum_header) + '\n')

    if not os.path.isdir(result_folder + '/analyzing/'):
        os.mkdir(result_folder + '/analyzing/')

    if not os.path.isdir(accumulated_reads):
        os.mkdir(accumulated_reads)

    if os.path.isdir(fastq_path):
        # print("yes")
        # print(os.listdir(fastq_path))
        folder_name_list = [name for name in os.listdir(fastq_path) if os.path.isdir(args.path + '/fastq_pass/' + name)]
        # print(folder_name_list)
        if len(folder_name_list) == 0:
            return ['no_subfolder_in_fastq_pass', 0]

        else:
            fastq_file_all = glob.glob(str(fastq_path + '/*/*.fastq'))
            # print(fastq_file_all)
            if len(fastq_file_all) >= 1:
                if not os.path.isdir(analyzing):
                    os.mkdir(analyzing)

                fastq_file_count = 0
                for name in folder_name_list:
                    fastq_list = glob.glob(str(fastq_path + '/' + name + '/*.fastq'))
                    if len(fastq_list) >= 1:
                        fastq_file_count += len(fastq_list)
                        mv_dir = analyzing + '/' + name + '/'
                        os.mkdir(mv_dir)
                        for file in fastq_list:
                            shutil.move(file, mv_dir)
                return ['find_new_fastq', fastq_file_count]
            else:
                return ['no_new_fastq_file', 0]
    else:
        return ['no_fastq_pass_folder', 0]


def bash_analyze(sample_path, fq_path, one_sample, thread, refer_seq, accumulated_reads, target_region, summary,
                 pooled_result):
    awk = """awk 'BEGIN {sum=0} {sum+=$3} END {print "mapped_depth", sum}' """
    awk2 = """awk 'NR%4==2{c++; l+=length($0)} END{print "read_number", c; print "base_number", l}'"""
    awk3 = """awk '{print $2, $1}' """
    awk4 = """awk '{if($2==0)print $1, "1"; else print}' """

    cmd = """cat {fq_path}/*fastq > {save}/result/{name}.fastq
             minimap2 -t {thread} -a -x map-ont -Y --MD {refer} {save}/result/{name}.fastq > {save}/result/{name}.sam \
                                    2>> {save}/result/{name}_alignment_summary.log
             samtools sort -@ {thread} {save}/result/{name}.sam -o {save}/result/{name}.bam >> {save}/result/{name}_alignment_summary.log 2>&1
             samtools index -@ {thread} {save}/result/{name}.bam >> {save}/result/{name}_alignment_summary.log
             samtools depth {save}/result/{name}.bam > {save}/result/{name}.coverage.txt
             cat {save}/result/{name}.fastq | {awk2} >> {save}/result/{name}_alignment_summary.log
             cat {save}/result/{name}.coverage.txt | {awk} >> {save}/result/{name}_alignment_summary.log
             grep -v "^@" {save}/result/{name}.sam | cut -f 3 | sort | uniq -c | sort -n | {awk3} >> {save}/result/{name}_alignment_summary.log
             cat {save}/result/{name}.fastq >> {accumulated_reads}/{name}.fastq
             cat {save}/result/{name}_alignment_summary.log | tail -n 5 | {awk4} > {save}/result/{name}_alignment_summary.tail5.log
             """.format(save=sample_path,
                        fq_path=fq_path,
                        name=one_sample,
                        thread=thread,
                        refer=refer_seq,
                        accumulated_reads=accumulated_reads,
                        awk=awk,
                        awk2=awk2,
                        awk3=awk3, awk4=awk4)
    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if a.stdout != b'':
        logging.info(a.stdout.decode('utf-8'))

    log_file = sample_path + '/result/' + one_sample + '_alignment_summary.tail5.log'
    coverage_file = sample_path + '/result/' + one_sample + '.coverage.txt'
    # count on target base
    on_target_base(coverage_file, target_region, log_file)
    print_result(one_sample, log_file, summary, pooled_result, target_region)


def individual_analysis(args, result_folder):
    main_analyzing = result_folder + '/analyzing/'
    analyzed = result_folder + '/analyzed_achieve'
    accumulated_reads = result_folder + '/analyzed_achieve/accumulated_reads'

    ctime = [name for name in os.listdir(main_analyzing) if os.path.isdir(main_analyzing + name)]
    if len(ctime) == 1:
        analyzing = main_analyzing + ctime[0] + '/'
    else:
        sys.stderr.write("ERROR!  No folder or more than 1 folder in %s\n" % main_analyzing)
        sys.exit(1)

    summary = main_analyzing + ctime[0] + '/' + ctime[0] + '_result.txt'
    pooled_result = result_folder + '/analyzed_achieve/pooled_result.txt'
    updated_result = result_folder + '/' + ctime[0] + '_result.txt'

    all_sample = [name for name in os.listdir(analyzing) if os.path.isdir(analyzing + name)]

    for one_sample in all_sample:
        sample_path = analyzing + one_sample
        os.mkdir(sample_path + '/result')

        if args.guppy_barcoder is not None:
            cmd = """mkdir {save}/fastq
                     mv {save}/*.fastq {save}/fastq
                     {gp} --require_barcodes_both_ends -i {save}/fastq -s {save} --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" \
                        -t {thread} >> {save}/result/{name}_alignment_summary.log
                     """.format(save=sample_path,
                                gp=args.guppy_barcoder,
                                name=one_sample,
                                thread=args.thread)
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

            gp_fastq_dir = sample_path + '/' + one_sample
            if os.path.isdir(gp_fastq_dir):
                fq_path = gp_fastq_dir
                bash_analyze(sample_path, fq_path, one_sample, args.thread, args.refer_seq, accumulated_reads,
                             args.target_region, summary, pooled_result)

        else:
            fq_path = sample_path
            bash_analyze(sample_path, fq_path, one_sample, args.thread, args.refer_seq, accumulated_reads,
                         args.target_region, summary, pooled_result)

    old_result = result_folder + '/*_result.txt'
    old_result = glob.glob(old_result)
    if len(old_result) >= 1:
        for file in old_result:
            shutil.move(file, analyzed)

    cmd = """mv {analyzing} {analyzed}""".format(analyzed=analyzed, analyzing=analyzing)
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    result_pool(pooled_result, updated_result, args.target_region)


def on_target_base(depth_file, target_region, log_file):
    # depth_file = './barcode13.coverage.txt'
    # target_region = '11049-11236,14306-14500,23124-23431,28086-28753'
    # log_file = './barcode13_alignment_summary.tail5.log'

    target_region = target_region.split(',')

    if len(target_region) == 0:
        logging.error("ERROR!  No --target_region provided\n"
                      "ERROR!  Please provide --target_region like: '11049-11236,14306-14500' \n")
        sys.exit(1)

    with open(depth_file, 'r') as infile, open(log_file, 'a') as outfile:
        depth_dir = {}
        for line in infile:
            n, p, c = line.strip().split()
            depth_dir[p.strip()] = c.strip()

        for region in target_region:
            region_split = re.split(r'_|-', region)
            region_l = int(region_split[0])
            region_r = int(region_split[1])

            region_base = 0

            for key in range(region_l, region_r+1):
                if str(key) in depth_dir:
                    region_base += int(depth_dir[str(key)])

            # if region_base == 0:
            #     region_base = 1
            outfile.writelines(str(region) + ' ' + str(region_base) + '\n')


def print_result(name, log_file, summary, pooled_result, target_region):
    target_region = target_region.split(',')
    with open(log_file, 'r') as infile, open(summary, 'a') as outfile, open(pooled_result, 'a') as outfile2:
        result_dir = {}
        for line in infile:
            key, value = line.strip().split()
            result_dir[key.strip()] = value.strip()

        barcode = name
        read_number = result_dir['read_number']
        base_number = result_dir['base_number']
        mapped_depth = result_dir['mapped_depth']
        mapped_record = result_dir['SARS-CoV-2']  # extract name from reference sequence
        unmapped_record = result_dir['*']

        ontarget_base = 0
        region_list = []
        for region in target_region:
            ontarget_base += int(result_dir[str(region)])
            region_list = region_list + [str(result_dir[str(region)])]

        total_record_number = int(mapped_record) + int(unmapped_record)
        percentage_of_mapped_record = round(int(mapped_record) / int(total_record_number), 4)
        percentage_of_mapped_base = round(int(mapped_depth) / int(base_number), 4)
        percentage_of_ontarget_base = round(int(ontarget_base) / int(mapped_depth), 4)

        result = [barcode, str(percentage_of_mapped_record), str(percentage_of_mapped_base),
                  str(percentage_of_ontarget_base), str(read_number), str(base_number), str(mapped_record),
                  str(unmapped_record), str(mapped_depth), str(ontarget_base)]

        result = result + region_list

        outfile.writelines('\t'.join(result) + '\n')
        outfile2.writelines('\t'.join(result) + '\n')


def put_back(args):
    input_path = args.path
    if args.save_path is None:
        result_folder = input_path + '/rtnano_result'
    else:
        result_folder = args.save_path

    if args.guppy_barcoder is not None:
        fastq_in = '/fastq/'
    else:
        fastq_in = '/'

    fastq_path = args.path + '/fastq_pass/'

    analyzed_folder = result_folder + '/analyzed_achieve/'

    if os.path.isdir(analyzed_folder):
        fastq_count = 0
        time_folder_list = [name for name in os.listdir(analyzed_folder) if os.path.isdir(analyzed_folder + name)]
        for time_folder in time_folder_list:
            if time_folder.startswith("20"):
                time_folder_dir = analyzed_folder + time_folder + '/'
                sub_folder_list = [name for name in os.listdir(time_folder_dir) if os.path.isdir(time_folder_dir + name)]
                for barcode_folder in sub_folder_list:
                    barcode_folder_dir = time_folder_dir + barcode_folder
                    fastq_dir = barcode_folder_dir + fastq_in

                    fastq_list = glob.glob(str(fastq_dir + '*.fastq'))
                    if len(fastq_list) >= 1:
                        mv_dir = fastq_path + barcode_folder + '/'
                        for file in fastq_list:
                            fastq_count += 1
                            shutil.move(file, mv_dir)
        print("%s    Finished! total %s fastq file are returned to <%s>\n" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                                                                              str(fastq_count), fastq_path))

    else:
        sys.stderr.write("ERROR!  analyzed_folder not exist <%s>\n"
                         "ERROR!  Please use the same command line as you generate the result, "
                         "together with --put_back\n" % analyzed_folder)
        sys.exit(1)


def result_pool(pooled_result, updated_result, target_region):
    target_region = target_region.split(',')
    df = pd.read_table(pooled_result)
    column_names = ["#barcode", "%_mapped_record", "%_mapped_base", "%_ontarget_base", "read_number", "total_base",
                    "mapped_record", "unmapped_record", "mapped_base", "ontarget_base"]

    column_names = column_names + target_region
    new_df = pd.DataFrame(columns=column_names)

    barcode_list = df['#barcode'].unique()

    for bc in barcode_list:
        bc_df = df[df['#barcode'].str.match(bc)]
        bc_df_sum = bc_df.sum(axis=0)

        if bc == 'unclassified':
            bc = 'unclassified99'

        bc_df_sum['#barcode'] = bc
        new_df.loc[-1] = bc_df_sum
        new_df.index = new_df.index + 1
        new_df = new_df.sort_index()

    change_dtype_lst = column_names[4:]
    for column in change_dtype_lst:
        new_df[column] = new_df[column].astype(str).astype(int)

    new_df['%_mapped_record'] = (new_df['mapped_record'] / (new_df['mapped_record'] + new_df['unmapped_record'])).round(4)
    new_df['%_mapped_base'] = (new_df['mapped_base'] / new_df['total_base']).round(4)
    new_df['%_ontarget_base'] = (new_df['ontarget_base'] / new_df['mapped_base']).round(4)
    new_df['sort'] = new_df['#barcode'].str.extract('(\d+)', expand=False).astype(int)
    new_df.sort_values('sort', inplace=True)
    new_df = new_df.drop('sort', axis=1)

    new_df.to_csv(updated_result, header=True, index=None, sep='\t', mode='w')
    with open(updated_result, 'r') as infile:
        for line in infile:
            logging.info(line.strip())


def cycle_run(args, result_folder, cycle_time):
    for cycle in range(cycle_time):
        result = get_fastq_file(args, result_folder)

        if result[0] == 'no_fastq_pass_folder':
            logging.info("%s    WARNING!    No fastq_pass folder in input folder, will check again in %s min" %
                         (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), args.interval_time))
            time.sleep(args.interval_time * 60)

        elif result[0] == 'no_subfolder_in_fastq_pass':
            logging.info("%s    No demultiplexed barcode folder in fastq_pass folder, will check again in %s min" %
                         (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), args.interval_time))
            time.sleep(args.interval_time * 60)

        elif result[0] == 'no_new_fastq_file':
            logging.info("%s    No new fastq file in fastq_pass folder, will check again in %s min" %
                         (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), args.interval_time))
            time.sleep(args.interval_time * 60)

        elif result[0] == 'find_new_fastq':
            logging.info("%s    -->Start analyzing %s new fastq file ...\n" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                                                                               str(result[1])))
            start = timer()

            individual_analysis(args, result_folder)

            end = timer()
            used_time = round(end - start)
            used_time_min = round(used_time / 60)

            if int(args.interval_time) > int(used_time_min):
                left_time = int(args.interval_time) - int(used_time_min)
                logging.info("\n%s    <--Analysis finished in %s s, next run in %s min" %
                             (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(used_time), str(left_time)))
                time.sleep(left_time * 60)
            else:
                logging.info("\n%s    <--Analysis finished in %s s, next run start now" %
                             (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(used_time)))

        else:
            logging.error("ERROR!  get_fastq_file return unrecognised value: %s\n"
                          "ERROR!  check get_fastq_file(args, result_folder)\n" % result)
            sys.exit(1)


def main():
    args = get_argparse()

    if args.put_back is True:
        print("\n%s    Program start ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        put_back(args)

    elif args.call_variant is True:
        print("\n%s    Program start ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print("--> Read fastq file from %s/analyzed_achieve/accumulated_reads/*.fastq" % args.save_path)

        if args.save_path is None:
            args.save_path = args.path + '/rtnano_result'

        fastq_regex = args.save_path + '/analyzed_achieve/accumulated_reads/*.fastq'

        fastq_file_all = glob.glob(str(fastq_regex))
        if len(fastq_file_all) >= 1:
            print("--> Detect %s fastq file" % len(fastq_file_all))
            print("--> Result saved in %s/snv/%s" % (args.save_path,
                                                     str(datetime.now().strftime("%Y%m%d_%H.%M.%S"))))

            save_path = args.save_path + '/snv/' + str(datetime.now().strftime("%Y%m%d_%H.%M.%S")) + '/'
            if not os.path.isdir(args.save_path + '/snv/'):
                os.mkdir(args.save_path + '/snv/')

            if not os.path.isdir(args.save_path + '/snv/' + str(datetime.now().strftime("%Y%m%d_%H.%M.%S")) + '/'):
                os.mkdir(args.save_path + '/snv/' + str(datetime.now().strftime("%Y%m%d_%H.%M.%S")) + '/')

            p = multiprocessing.Pool(args.thread)
            for file in fastq_file_all:
                p.apply_async(call_variant.call, args=(file, save_path, args.refer_seq, args.alle_freq))

            p.close()
            p.join()
            print("%s    Variant calling finished!\n" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        else:
            sys.stderr.write("ERROR! No fastq file detected in %s/analyzed_achieve/accumulated_reads/\n" % args.save_path)
            sys.exit(0)

    else:
        result_folder = prepare_env(args)

        ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
        log_file = result_folder + '/' + ctime + '_rt_nano.log'
        logging.basicConfig(level=logging.DEBUG,
                            format='%(message)s',
                            filename=log_file,
                            filemode='w')
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
        logging.info("\n%s    Program start ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        logging.info("--> Working thread: %s" % args.thread)
        logging.info("--> Reference genome: %s" % args.refer_seq)
        logging.info("--> Target regions: %s" % args.target_region)
        logging.info("--> Result saved in %s" % result_folder)
        logging.info("--> Log saved in    %s" % log_file)

        cycle_time = round(args.run_time * 60 / args.interval_time)
        cycle_run(args, result_folder, cycle_time)
        logging.info("\n%s    Run finished! Result and log saved in %s" %
                     (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), result_folder))


if __name__ == '__main__':
    main()
