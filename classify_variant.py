import logging
import multiprocessing


def get_snv(test, vcf_file="", rule_file=""):
    vcf_file = "/Users/bic/Desktop/codes/github/RTNano/RTNano/required_files/barcode18.snp.flt.vcf"
    rule_file = "/Users/bic/Desktop/codes/github/RTNano/RTNano/required_files/rules_test.txt"

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

        for variant in unique_variant_set:
            # print(variant)
            if variant in rule_file_mutations:
                # print(variant + "\t" + rule_file_dic[variant])
                detect_variant_per_strain_count[rule_file_dic[variant]] += 1
            else:
                logging.info("WARNING!  unknown variant type!\n" + variant + "\tunknown")

        # result = {}
        # for strain in set(rule_file_strain):
        #     result[str(strain)] = str(str(detect_variant_per_strain_count[strain]) + "/" + str(variant_number_per_strain[strain]))
        # print(result)

        result = []
        for strain in set(rule_file_strain):
            result.append("[" + str(strain) + "](" + str(
                str(detect_variant_per_strain_count[strain]) + "/" + str(variant_number_per_strain[strain])) + ")")

        # print("\t".join(result))
        return "\t".join(result)


# print(get_snv())

p = multiprocessing.Pool(5)
result = []

for i in range(10):
    res = p.apply_async(get_snv, args=(1,))
    result.append(res)

p.close()
p.join()

for res in result:
    print(res.get())
