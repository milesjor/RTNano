import logging
import subprocess


def get_identity(sample_path, fq_path, one_sample, thread, refer_seq, accumulated_reads, target_amplicon, summary,
                 pooled_result, alignment_identity, covered_region):
    awk = """awk '{if($6 == "'${ref}'")print}' | \
             awk 'BEGIN{a=0;b=0;c=0}{a+=$10/$11;b++;c+=($9-$8+1)/$7}END{print "'$ref'", a/b"/"c/b"/"b}' """

    awk2 = """awk 'NR%4==2{c++; l+=length($0)} END{print "read_number", c; print "base_number", l}'"""

    awk_filter = """awk '{if($6 == "'${ref}'")print}' | \
                    awk 'BEGIN{a=0;b=0;c=0;d=0;e=0}{a=$10/$11;b=($9-$8+1)/$7; \
                    if(a>="'${ident}'" && b>="'${cover}'"){d+=a;e+=b;c++}} \
                    END{if(c>=1) print "'$ref'", d/c"/"e/c"/"c}' """

    cmd = """cat {fq_path}/*fastq > {save}/result/{name}.fastq
             minimap2 -x map-ont -c -t {thread} {refer} {save}/result/{name}.fastq 1> {save}/result/{name}.paf \
                                    2>> {save}/result/{name}_alignment_summary.log
             refNames=$(cat {save}/result/{name}.paf | cut -f6 | sort | uniq)
             ident={ident}
             cover={cover}
             
             for ref in $refNames ; do
             # cat {save}/result/{name}.paf | {awk} >> {save}/result/{name}_alignment_summary.result.log
             cat {save}/result/{name}.paf | {awk_filter} >> {save}/result/{name}_alignment_summary.result.log
             done
             
             cat {save}/result/{name}.fastq | {awk2} >> {save}/result/{name}_alignment_summary.result.log
             cat {save}/result/{name}.fastq >> {accumulated_reads}/{name}.fastq
             """.format(save=sample_path,
                        fq_path=fq_path,
                        name=one_sample,
                        thread=thread,
                        refer=refer_seq,
                        accumulated_reads=accumulated_reads,
                        awk=awk,
                        awk2=awk2,
                        awk_filter=awk_filter,
                        ident=alignment_identity,
                        cover=covered_region)

    a = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if a.stdout != b'':
        logging.info(a.stdout.decode('utf-8'))

    log_file = sample_path + '/result/' + one_sample + '_alignment_summary.result.log'
    clean_result(one_sample, log_file, summary, pooled_result, target_amplicon)


def clean_result(name, log_file, summary, pooled_result, target_amplicon):

    with open(log_file, 'r') as infile, open(summary, 'a') as outfile, open(pooled_result, 'a') as outfile2:
        result_dir = {}
        for line in infile:
            key, value = line.strip().split()
            result_dir[key.strip()] = value.strip()

        barcode = name
        read_number = result_dir['read_number']
        base_number = result_dir['base_number']

        amplicon_record = []
        for name in target_amplicon:
            if name in result_dir:
                amplicon_record = amplicon_record + [str(result_dir[name])]
            else:
                amplicon_record = amplicon_record + ["0/0/0"]

        result = [barcode, str(read_number), str(base_number)] + amplicon_record

        outfile.writelines('\t'.join(result) + '\n')
        outfile2.writelines('\t'.join(result) + '\n')
