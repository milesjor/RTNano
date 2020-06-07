import re
import subprocess
import os
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
             bcftools mpileup --threads {thread} -Q 10 -Ou -A -f {refer} {save}/{name}.sv.sorted.bam 2>>{save}/{name}_alignment_summary.log | \
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
        print("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
    if a.stderr != b'':
        print("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))


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
            print("stdout from " + save_path + " :\n" + a.stdout.decode('utf-8'))
        if a.stderr != b'':
            print("stderr from " + save_path + " :\n" + a.stderr.decode('utf-8'))

    else:
        print("WARNING! VCF file do not exist -> %s" % vcf)


def call(fastq, save_path, refer, alle_freq, thread='1'):
    variant_calling(fastq, save_path, thread, refer, ax='map-ont')
    add_pcent(fastq, save_path, alle_freq)
    print("%s    %s finished!" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), fastq))
