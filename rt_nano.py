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
import logging
import multiprocessing
import pandas as pd
from _version import version
from datetime import datetime
from timeit import default_timer as timer
import call_variant
import identity_analysis
import summarize_result


def get_argparse():
    parser = argparse.ArgumentParser(description='Real-Time analysis of Nanopore data for Covid-19 sequencing.')
    parser.add_argument('-p', '--path', type=str, required=True, help='path/to/nanopore_result_folder')
    parser.add_argument('-s', '--save_path', type=str,
                        help='path/to/saved_folder Default: rtnano_result in -p PATH folder')
    parser.add_argument('-r', '--refer_seq', type=validate_file,
                        help='path/to/reference_genome.fa, default is using amplicon.fa in program folder')
    parser.add_argument('-t', '--thread', type=int, default='1', help='working thread [1]')
    parser.add_argument('-T', '--interval_time', type=int, default='1',
                        help='interval time for scanning minknow folder in minutes [1]')
    parser.add_argument('-g', '--guppy_barcoder', type=validate_file,
                        help='Optional: path/to/guppy_barcoder, when offering this parameter, it will do additional '
                             'demultiplexing using guppy_barcoder --require_barcodes_both_ends --trim_barcodes')
    parser.add_argument('-k', '--barcode_kits', type=str,
                        help='barcode kits used, e.g. "EXP-NBD114 EXP-NBD104" it is required when providing -g/--guppy_barcoder')
    parser.add_argument('--run_time', type=int, default='48',
                        help='total run time in hours [48]')
    parser.add_argument('--align_identity', type=float, default='0.89',
                        help='minimum [0.89] alignment identity when considering an effective alignment')
    parser.add_argument('--covered_pcent', type=float, default='0.96',
                        help='minimum covering [0.96] of amplicon length when considering an effective alignment')
    parser.add_argument('--housekeep_gene', type=str, default='ACTB_263bp',
                        help='the amplicon name of housekeeping gene in the sequencing, default: [ACTB_263bp]')
    parser.add_argument('--rule_file', type=validate_file,
                        help='the rule_file of strain variants, default is '
                             'using RTNano/required_files/variant_rules.txt in program folder')
    parser.add_argument('--ntc', type=str,
                        help='barcode number of no template control in the sequencing, e.g. barcode96')
    parser.add_argument('--resume', action='store_true',
                        help='resume the unexpectedly interrupted analysis. Please use the same [-p] [-s] as before')
    parser.add_argument('--put_back', action='store_true',
                        help='return fastq file to their original fastq_pass folder. Please use the same [-p] [-s] [-g]'
                             ' as you generated the result, together with --put_back')
    parser.add_argument('--rt_variant', action='store_true',
                        help='real-time call variants using samtools and filter by alleic frequency (>=0.5). ')
    parser.add_argument('--call_variant', action='store_true',
                        help='call variants using samtools and filter by alleic frequency (>=0.5). '
                             'Please use the same [-p] [-s] '
                             'as you generated the result. It uses fastq file in analyzed_achieve/accumulated_reads '
                             'folder. If you want to use it for your own data, '
                             'please put fastq file in this folder, one sample one fastq file')
    parser.add_argument('--alle_freq', type=float, default='0.5', help='filter SNVs by allelic frequency [0.5]')
    parser.add_argument('-v', '--version', action='version', version=version, help='show the current version')
    args = parser.parse_args()
    if args.refer_seq is None:
        args.refer_seq = sys.path[0] + '/amplicon.fa'
    if args.rule_file is None:
        args.rule_file = sys.path[0] + '/required_files/variant_rules.txt'
    if args.guppy_barcoder is not None:
        if args.barcode_kits is None:
            sys.stderr.write("ERROR!  Please provide -k/--barcode_kits\n")
            sys.exit(1)
        else:
            kit_list = barcode_kit_list()
            for kit in str(args.barcode_kits).split(" "):
                if kit not in kit_list:
                    print_kit_list = "\n".join(kit_list)
                    sys.stderr.write("ERROR!  -k/--barcode_kits is not in kit list:\n\n%s\n\n" % print_kit_list)
                    sys.exit(1)
    return args


def validate_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def barcode_kit_list():
    choices = ["EXP-NBD103", "EXP-NBD104", "EXP-NBD114", "EXP-NBD196", "EXP-PBC001", "EXP-PBC096",
               "OND-SQK-LP0096M", "OND-SQK-LP0096S", "OND-SQK-LP1152S", "OND-SQK-LP9216",
               "SQK-16S024", "SQK-LWB001", "SQK-PBK004", "SQK-PCB109", "SQK-RAB201", "SQK-RAB204",
               "SQK-RBK001", "SQK-RBK004", "SQK-RBK096", "SQK-RLB001", "SQK-RPB004",
               "VSK-VMK001", "VSK-VMK002"]
    return choices


def check_tool(tool):
    """Check whether `tool` is on PATH and marked as executable."""
    return shutil.which(tool) is not None


def prepare_env(args):
    input_path = args.path

    if args.save_path is None:
        result_folder = input_path + '/rtnano_result'
    else:
        result_folder = args.save_path

    if os.path.isdir(result_folder):
        if args.resume is not True:
            sys.stderr.write("WARNING!  Result folder exist! <%s> \n"
                             "WARNING!  To avoid overwriting in existing data, "
                             "Please use a new path and restart the program.\n"
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
        with open(pooled_result, 'w') as outfile, open(args.refer_seq) as infile:
            target_amplicon = []

            for line in infile:
                if line.startswith(">"):
                    ref_name = line.strip(">").split()[0]
                    target_amplicon.append(ref_name)

            header = ["#barcode", "read_number", "total_base"] + target_amplicon
            outfile.writelines('\t'.join(header) + '\n')

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
            cmd = """find {gz_fastq} -name "*.fastq.gz" -print0 | xargs -0 gunzip""".format(gz_fastq=fastq_path)
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

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


def individual_analysis(args, result_folder, target_amplicon):
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

            # cmd = """mkdir {save}/fastq
            #                      mv {save}/*.fastq {save}/fastq
            #                      {gp} --require_barcodes_both_ends -i {save}/fastq -s {save} \
            #                         --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg" \
            #                         -t {thread} --trim_barcodes >> {save}/result/{name}_alignment_summary.log

            cmd = """mkdir {save}/fastq
                     mv {save}/*.fastq {save}/fastq
                     {gp} --require_barcodes_both_ends -i {save}/fastq -s {save} \
                        --barcode_kits "{barcode_kits}" \
                        -t {thread} --trim_barcodes >> {save}/result/{name}_alignment_summary.log
                     """.format(save=sample_path,
                                gp=args.guppy_barcoder,
                                name=one_sample,
                                thread=args.thread,
                                barcode_kits=args.barcode_kits)
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

            gp_fastq_dir = sample_path + '/' + one_sample
            if os.path.isdir(gp_fastq_dir):
                fq_path = gp_fastq_dir
                identity_analysis.get_identity(sample_path, fq_path, one_sample, args.thread, args.refer_seq,
                                               accumulated_reads, target_amplicon, summary, pooled_result,
                                               args.align_identity, args.covered_pcent)

        else:
            fq_path = sample_path
            identity_analysis.get_identity(sample_path, fq_path, one_sample, args.thread, args.refer_seq,
                                           accumulated_reads, target_amplicon, summary, pooled_result,
                                           args.align_identity, args.covered_pcent)

    old_result = result_folder + '/*_result.txt'
    old_result = glob.glob(old_result)
    if len(old_result) >= 1:
        for file in old_result:
            shutil.move(file, analyzed)

    cmd = """mv {analyzing} {analyzed}""".format(analyzed=analyzed, analyzing=analyzing)
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

    result_pool(pooled_result, updated_result, target_amplicon, args.housekeep_gene, args.ntc)


def result_pool(pooled_result, updated_result, target_amplicon, housekeeping, NTC):
    column_names = ["#barcode", "read_number", "total_base"] + target_amplicon

    df = pd.read_table(pooled_result)
    new_df = pd.DataFrame(columns=column_names)

    barcode_list = df['#barcode'].unique()

    for bc in barcode_list:
        bc_df = df[df['#barcode'].str.match(bc)]
        bc_df_sum = bc_df.sum(axis=0)

        for amplicon in target_amplicon:
            part_df = bc_df[amplicon]
            part_df = part_df.str.split('/', expand=True).astype(float)

            part_df[0] = part_df[0] * part_df[2]
            part_df[1] = part_df[1] * part_df[2]
            part_df_sum = part_df.sum(axis=0)

            if part_df_sum[2] == 0:
                bc_df_sum[amplicon] = 'N/A'
            else:
                part_df_sum[0] = round(part_df_sum[0] / part_df_sum[2] * 100, 1)
                part_df_sum[1] = round(part_df_sum[1] / part_df_sum[2] * 100, 1)

                bc_df_sum[amplicon] = part_df_sum[0].astype(str) + '/' + part_df_sum[1].astype(str) + '/' + \
                                      part_df_sum[2].astype(int).astype(str)

        if bc == 'unclassified':
            bc = 'unclassified99'

        bc_df_sum['#barcode'] = bc
        new_df.loc[-1] = bc_df_sum
        new_df.index = new_df.index + 1
        new_df = new_df.sort_index()

    new_df['sort'] = new_df['#barcode'].str.extract('(\d+)', expand=False).astype(int)
    new_df.sort_values('sort', inplace=True)
    new_df = new_df.drop('sort', axis=1)

    new_df.to_csv(updated_result, header=True, index=None, sep='\t', mode='w')
    with open(updated_result, 'r') as infile:
        for line in infile:
            logging.info(line.strip())

    # add result label

    logging.info("\n- Clean result -\n")
    summarize_result.summarize(updated_result, target_amplicon, housekeeping, NTC)


def cycle_run(args, result_folder, cycle_time, target_amplicon):
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
            logging.info("%s    -->Start analyzing %s new fastq file ...\n" %
                         (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(result[1])))
            start = timer()

            individual_analysis(args, result_folder, target_amplicon)
            # ---------------------------------
            if args.rt_variant is True:
                variant_call(args)
            # ---------------------------------
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
    analyzing_folder = result_folder + '/analyzing/'

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

        if os.path.isdir(analyzing_folder):
            folder_list = [name for name in os.listdir(analyzing_folder) if os.path.isdir(analyzing_folder + name)]
            if len(folder_list) >= 1:
                for time_folder in folder_list:
                    if time_folder.startswith("20"):
                        time_folder_dir = analyzing_folder + time_folder + '/'
                        sub_folder_list = [name for name in os.listdir(time_folder_dir) if
                                           os.path.isdir(time_folder_dir + name)]
                        for barcode_folder in sub_folder_list:
                            barcode_folder_dir = time_folder_dir + barcode_folder

                            if os.path.isdir(barcode_folder_dir + '/fastq'):
                                fastq_in = '/fastq/'
                            else:
                                fastq_in = '/'

                            fastq_dir = barcode_folder_dir + fastq_in
                            fastq_list = glob.glob(str(fastq_dir + '*.fastq'))
                            if len(fastq_list) >= 1:
                                mv_dir = fastq_path + barcode_folder + '/'
                                for file in fastq_list:
                                    fastq_count += 1
                                    shutil.move(file, mv_dir)

        print("%s    Finished! total %s fastq file are returned to <%s>\n" %
              (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(fastq_count), fastq_path))

    else:
        sys.stderr.write("ERROR!  analyzed_folder not exist <%s>\n"
                         "ERROR!  Please use the same command line as you generate the result, "
                         "together with --put_back\n" % analyzed_folder)
        sys.exit(1)


def variant_call(args):
    if args.save_path is None:
        args.save_path = args.path + '/rtnano_result'

    fastq_regex = args.save_path + '/analyzed_achieve/accumulated_reads/*.fastq'

    fastq_file_all = glob.glob(str(fastq_regex))
    if len(fastq_file_all) >= 1:
        logging.info("\n--> Variant calling detect %s fastq file" % len(fastq_file_all))
        logging.info("--> Variant result saved in %s/snv/%s\n" % (args.save_path,
                                                 str(datetime.now().strftime("%Y%m%d_%H.%M.%S"))))

        save_path = args.save_path + '/snv/' + str(datetime.now().strftime("%Y%m%d_%H.%M.%S")) + '/'
        if not os.path.isdir(args.save_path + '/snv/'):
            os.mkdir(args.save_path + '/snv/')

        if not os.path.isdir(args.save_path + '/snv/' + str(datetime.now().strftime("%Y%m%d_%H.%M.%S")) + '/'):
            os.mkdir(args.save_path + '/snv/' + str(datetime.now().strftime("%Y%m%d_%H.%M.%S")) + '/')

        p = multiprocessing.Pool(args.thread)

        for file in fastq_file_all:
            p.apply_async(call_variant.call, args=(file, save_path, args.refer_seq, args.alle_freq, args.rule_file))

        p.close()
        p.join()

        save_log_file = save_path + "/variant_summary.txt"
        with open(save_log_file, "r") as save_file:
            for line in sorted(save_file):
                logging.info(line.strip())
        logging.info("\n%s    Variant calling finished!" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def main():
    args = get_argparse()
    if check_tool("minimap2") is not True:
        sys.exit("ERROR! Executable minimap2 is not found!\n"
                 "ERROR! Please install minimap2 -> `conda install minimap2=2.11`")

    if args.put_back is True:
        print("\n%s    Program start ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        put_back(args)

    elif args.call_variant is True:
        if check_tool("samtools") is not True:
            sys.exit("ERROR! Executable samtools is not found!\n"
                     "ERROR! Please install samtools -> `conda install samtools=1.9`")
        if check_tool("bcftools") is not True:
            sys.exit("ERROR! Executable bcftools is not found!\n"
                     "ERROR! Please install bcftools -> `conda install bcftools=1.9`")

        if args.save_path is None:
            args.save_path = args.path + '/rtnano_result'

        result_folder = args.save_path
        ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
        log_file = result_folder + '/' + ctime + '_rt_nano.log'
        logging.basicConfig(level=logging.DEBUG,
                            format='%(message)s',
                            filename=log_file,
                            filemode='w')
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

        fastq_regex = args.save_path + '/analyzed_achieve/accumulated_reads/*.fastq'
        fastq_file_all = glob.glob(str(fastq_regex))

        if len(fastq_file_all) >= 1:
            logging.info("\n%s    Program start ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            logging.info("--> Read fastq file from %s/analyzed_achieve/accumulated_reads/*.fastq" % args.save_path)
            variant_call(args)

        else:
            sys.stderr.write("ERROR! No fastq file detected in %s/analyzed_achieve/accumulated_reads/\n" %
                             args.save_path)
            sys.exit(0)

    else:
        if args.rt_variant is True:
            if check_tool("samtools") is not True:
                sys.exit("ERROR! Executable samtools is not found!\n"
                         "ERROR! Please install samtools -> `conda install samtools=1.9`")
            if check_tool("bcftools") is not True:
                sys.exit("ERROR! Executable bcftools is not found!\n"
                         "ERROR! Please install bcftools -> `conda install bcftools=1.9`")

        result_folder = prepare_env(args)

        with open(args.refer_seq) as infile:
            target_amplicon = []

            for line in infile:
                if line.startswith(">"):
                    ref_name = line.strip(">").split()[0]
                    target_amplicon.append(ref_name)

        ctime = datetime.now().strftime("%Y%m%d_%H.%M.%S")
        log_file = result_folder + '/' + ctime + '_rt_nano.log'
        logging.basicConfig(level=logging.DEBUG,
                            format='%(message)s',
                            filename=log_file,
                            filemode='w')
        logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
        logging.info("\n%s    Program start ..." % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        logging.info("--> Working thread:  %s" % args.thread)
        logging.info("--> Reference genome: %s" % args.refer_seq)
        logging.info("--> Target amplicon: %s" % target_amplicon)
        logging.info("--> Result saved in:  %s" % result_folder)
        logging.info("--> Log saved in:     %s" % log_file)
        logging.info("--> housekeep gene:   %s" % args.housekeep_gene)
        logging.info("--> No template control:        %s" % args.ntc)
        logging.info("--> Covered region filter:      %s" % args.covered_pcent)
        logging.info("--> Alignment identity filter:  %s" % args.align_identity)
        logging.info("--> Positive sample confidence level: POS_3 > POS_2 > POS_1 > POS_0")

        cycle_time = round(args.run_time * 60 / args.interval_time)
        cycle_run(args, result_folder, cycle_time, target_amplicon)
        logging.info("\n%s    Run finished! Result and log saved in %s" %
                     (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), result_folder))


if __name__ == '__main__':
    main()
