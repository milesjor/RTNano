import re
import subprocess
import os
import logging
from datetime import datetime


def variant_calling(fastq, save_path, thread, refer, ax='map-ont'):

    file_name = re.split(r'/', fastq)[-1]
    file_name = re.split(r'.fastq', file_name)[0]
    save_path = save_path + '/' + file_name

    count_number = """awk 'BEGIN{n=0}{if($3>=3)n+=1}END{print n}' """
    cmd = """# alignment
             mkdir {save}
             mkdir {save}/vcf
             minimap2 -t {thread} -a -x {ax} -Y --MD {refer} {fastq} > {save}/{name}.sv.sam 2>>{save}/{name}_alignment_summary.log
             echo "====== Alignment summary ======" >> {save}/{name}_alignment_summary.log
             grep -v "^@" {save}/{name}.sv.sam | cut -f 3 | sort | uniq -c | sort -n >> {save}/{name}_alignment_summary.log
             echo "===============================" >> {save}/{name}_alignment_summary.log
             samtools sort -@ {thread} {save}/{name}.sv.sam -o {save}/{name}.sv.sorted.bam >> {save}/{name}_alignment_summary.log
             samtools index -@ {thread} {save}/{name}.sv.sorted.bam
             # SNP calling
             bcftools mpileup --threads {thread} -d 500 -Q 10 -Ou -A -f {refer} {save}/{name}.sv.sorted.bam 2>>{save}/{name}_alignment_summary.log | \
             bcftools call --threads {thread} --ploidy 1 -Ou -mv 2>>{save}/{name}_alignment_summary.log | \
             bcftools norm --threads {thread} -Ou -f {refer} 2>>{save}/{name}_alignment_summary.log | \
             bcftools annotate --set-id {name} -Ov 1> {save}/vcf/{name}.snp.vcf
             # bcftools filter -s LowQual -e '%QUAL<10 || INFO/DP<5 || INFO/IMF<0.5' > {save}/vcf/{name}.snp.vcf
             samtools depth {save}/{name}.sv.sorted.bam > {save}/{name}.coverage.txt
             cat {save}/{name}.coverage.txt | {count_number} | xargs echo {name} > {save}/{name}.coverage.3plus.txt
             """.format(save=save_path,
                        thread=thread,
                        ax=ax,
                        refer=refer,
                        fastq=fastq,
                        name=file_name,
                        count_number=count_number)

    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if a.stdout != b'':
        logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))


def add_pcent(fastq, save_path, alle_freq):

    file_name = re.split(r'/', fastq)[-1]
    file_name = re.split(r'.fastq', file_name)[0]
    save_path = save_path + '/' + file_name

    vcf = save_path + '/vcf/' + file_name + '.snp.vcf'
    outvcf = save_path + '/vcf/' + file_name + '.snp.flt.vcf.tmp'
    flt_vcf = save_path + '/vcf/' + file_name + '.snp.flt.vcf'

    if os.path.exists(vcf):
        with open(vcf, 'r') as infile, open(outvcf, 'w') as outfile1:
            for line in infile:
                if line.startswith('##'):
                    outfile1.writelines(line)
                elif line.startswith('#'):
                    outfile1.writelines(
                        """##INFO=<ID=VARP,Number=1,Type=Float,Description="Variant Percentage: Variant allele number divided by total allele number">\n""")
                    outfile1.writelines(
                        """##INFO=<ID=ALTC,Number=1,Type=Integer,Description="alt-read count from DP4">\n""")
                    outfile1.writelines(
                        """##INFO=<ID=REFC,Number=1,Type=Integer,Description="ref-read count from DP4">\n""")
                    outfile1.writelines(line)

                else:
                    line = re.split(r'\t', line)
                    info = line[7]
                    dp4 = re.split(r'DP4=|;MQ=', info)[1]
                    dp4 = re.split(r',', dp4)
                    var_reads = int(dp4[2]) + int(dp4[3])
                    ref_reads = int(dp4[0]) + int(dp4[1])
                    total_coverage = var_reads + ref_reads
                    var_percent = round(var_reads / total_coverage, 2)
                    new_info = info + ';ALTC=' + str(var_reads) + ';REFC=' + str(ref_reads) + ';VARP=' + str(
                        var_percent)
                    revised_line = tuple(line[0:7] + [new_info] + line[8:])
                    outfile1.writelines('\t'.join(revised_line))

        cmd = """cat {outvcf} | bcftools filter -s LowQual -e \
                    '%QUAL<10 || INFO/DP<3 || INFO/IMF<0.5 || INFO/VARP<{alle_freq} || INFO/ALTC<3' > {flt_vcf}
                 rm {outvcf}
                 """.format(outvcf=outvcf,
                            flt_vcf=flt_vcf,
                            alle_freq=alle_freq)
        a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if a.stdout != b'':
            logging.info("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
        if a.stderr != b'':
            logging.info("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    else:
        logging.info("WARNING! VCF file do not exist -> %s" % vcf)


def get_snv(vcf_file, rule_file, save_path, file_name):
    # vcf_file = "/Users/bic/Desktop/codes/github/RTNano/RTNano/required_files/barcode18.snp.flt.vcf"
    # rule_file = "/Users/bic/Desktop/codes/github/RTNano/RTNano/required_files/rules_test.txt"

    with open(vcf_file, "r") as vcf_file:
        variant_list = []  # list()
        for line in vcf_file:
            if not line.startswith("#"):
                # print(line.split("\t"))
                if not line.split("\t")[7].startswith("INDEL"):
                    variant = line.split("\t")[0] + ":" + line.split("\t")[1] + ":" + line.split("\t")[3].upper() \
                              + ">" + line.split("\t")[4].upper()
                    variant_list.append(variant)
                else:
                    # print(line)
                    ref_seq = line.split("\t")[3].upper()
                    call_seq = line.split("\t")[4].upper()
                    # print(ref_seq[0])
                    # print(call_seq[0])
                    if len(ref_seq) == 1 or len(call_seq) == 1:
                        if ref_seq[0] == call_seq[0]:
                            if len(ref_seq) > len(call_seq):
                                deletion_length = int(len(ref_seq)) - int(len(call_seq))
                                if deletion_length == 1:
                                    variant = line.split("\t")[0] + ":" + str(int(line.split("\t")[1]) + 1) + ":del" + str(ref_seq[-1])
                                    variant_list.append(variant)

                                else:
                                    variant = line.split("\t")[0] + ":" + str(int(line.split("\t")[1]) + 1) + "_" + str(
                                        int(line.split("\t")[1]) + deletion_length) + ":del"
                                    variant_list.append(variant)

                            elif len(ref_seq) < len(call_seq):
                                variant = line.split("\t")[0] + ":" + str(int(line.split("\t")[1])) + "_" + str(
                                    int(line.split("\t")[1]) + 1) + ":ins" + str(call_seq[1:])
                                variant_list.append(variant)

                        else:
                            logging.info("\tWARNING!  unknown variant type!")
                            logging.info(line)
                    else:
                        logging.info("\tWARNING!  unknown variant type!")
                        logging.info(line)

        unique_variant_set = set(variant_list)
        # print(unique_variant_set)

        rule_file_dic = {}
        rule_file_mutations = set()
        rule_file_strain = list()
        with open(rule_file, "r") as rule_file:
            for line in rule_file:
                # print(line.strip().split('\t', 1))
                strain, mutation = line.strip().split('\t', 1)
                rule_file_dic[mutation] = strain.strip()
                rule_file_mutations.add(mutation)
                rule_file_strain.append(strain)
        # print(rule_file_dic)
        # print(rule_file_mutations)

        variant_number_per_strain = {}
        for element in set(rule_file_strain):
            variant_number_per_strain[element] = rule_file_strain.count(element)
        # print(variant_number_per_strain)

        detect_variant_per_strain_count = {}
        for strain in set(rule_file_strain):
            detect_variant_per_strain_count[strain] = 0
        # print("aaaaaaaaa")
        # print(detect_variant_per_strain_count)

        save_unique_snv_file = save_path + "/all_unique_variant.txt"
        with open(save_unique_snv_file, "a") as save_file:
            # save_file.writelines("test")
            for variant in unique_variant_set:
                # print(variant)
                save_file.writelines(file_name + "\t" + variant + "\n")

                if variant in rule_file_mutations:
                    # print(variant + "\t" + rule_file_dic[variant])
                    detect_variant_per_strain_count[rule_file_dic[variant]] += 1
                # else:
                #     logging.info("WARNING!  unknown variant type!\n" + variant + "\tunknown")

        # result = {}
        # for strain in set(rule_file_strain):
        #     result[str(strain)] = str(str(detect_variant_per_strain_count[strain]) + "/" + str(variant_number_per_strain[strain]))

        result = []
        for strain in set(rule_file_strain):
            result.append("[" + str(strain) + "](" + str(
                str(detect_variant_per_strain_count[strain]) + "/" + str(variant_number_per_strain[strain])) + ")")

        return "\t".join(result)


def call(fastq, save_path, refer, alle_freq, rule_file, thread='1'):
    variant_calling(fastq, save_path, thread, refer, ax='map-ont')
    add_pcent(fastq, save_path, alle_freq)

    file_name = re.split(r'/', fastq)[-1]
    file_name = re.split(r'.fastq', file_name)[0]
    save_path_2 = save_path + '/' + file_name
    vcf_file = save_path_2 + '/vcf/' + file_name + '.snp.vcf'

    variant_classify = get_snv(vcf_file, rule_file, save_path, file_name)
    result = str(file_name + "\t" + variant_classify)

    save_log_file = save_path + "/variant_summary.txt"
    with open(save_log_file, "a") as save_file:
        # save_file.writelines("test")
        save_file.writelines(result + "\n")
    # logging.info("%s    %s finished!" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), fastq))
