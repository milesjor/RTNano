import pandas as pd
import logging


def rule():
    rules = {
        'POS_3': '3 region_alignment_record >= 50',
        'POS_2': '2 region_alignment_record >= 50',
        'POS_1': '1 region_alignment_record >= 50'
                'OR 2 region_alignment_record >= 20'
                'OR 3 region_alignment_record >= 10'
                'OR 4 region_alignment_record >= 5',
        'POS_0': 'ACTB house_keeping_record >= 1000'
                'AND any region_alignment_record >= 1',
        'NEG': 'ACTB house_keeping_record >= 1000'
                'AND all region_alignment_record = 0',
        'UNK': 'ACTB house_keeping_record < 1000'
                'AND all region_alignment_record = 0',
    }
    return rules


def summarize(updated_result, target_amplicon, housekeeping, NTC):
    # updated_result = '/Users/bic/OneDrive/seq_data/rtnano_Gerardo/20220324_11.16.15_result.txt'
    # target_amplicon = ['amplicon_1', 'amplicon_2', 'amplicon_3', 'amplicon_4', 'amplicon_5', 'ACTB_263bp']
    # housekeeping = 'ACTB_263bp'
    # NTC = 'barcode99'

    # print("target amplicon", target_amplicon)
    target_amplicon_no_hp = [n for n in target_amplicon if n != housekeeping]
    # print("target amplicon without housekeeping", target_amplicon_no_hp)

    outfile_header = '\t'.join(['\n\n#barcode', 'Result', 'read', 'base'] + target_amplicon) + '\n'
    # print(outfile_header)

    df = pd.read_table(updated_result)
    df = df.fillna('0/0/0')
    # print(df)

    summarized_result = [outfile_header]

    # get data of NTC
    ntc_record_dir = {}
    if NTC is not None:
        bc_df = df[df['#barcode'].str.match(NTC)]
        if bc_df.empty is True:
            logging.info("WARNING!   can not find NTC barcode data, ignore NTC in below result.")
            print("WARNING!   can not find NTC barcode data, ignore NTC in below result.")
            for amplicon in target_amplicon:
                ntc_record_dir[amplicon] = 0
        else:
            for amplicon in target_amplicon:
                part_df = bc_df[amplicon]
                part_df = part_df.str.split('/', expand=True).astype(float)
                # print(amplicon, "\n", part_df)
                ntc_record_dir[amplicon] = part_df[2].astype(int).values[0]
    else:
        for amplicon in target_amplicon:
            ntc_record_dir[amplicon] = 0

    barcode_list = df['#barcode'].unique()
    # barcode_list = ['barcode38']
    for bc in barcode_list:
        bc_df = df[df['#barcode'].str.match(bc)]
        # print(bc_df)
        record_dir = {}
        this_result_end = []

        for amplicon in target_amplicon:
            part_df = bc_df[amplicon]
            part_df = part_df.str.split('/', expand=True).astype(float)
            # print(amplicon, "\n", part_df)
            record_dir[amplicon] = part_df[2].astype(int).values[0]
            this_result_end.append(record_dir[amplicon])

        if housekeeping is None:
            housekeeping = "no_hp"
            record_dir[housekeeping] = 0

        # print("\nbarcode", bc)
        # print(record_dir)
        # print(this_result_end)

        cov2_50_plus = 0
        cov2_20_50 = 0
        cov2_10_20 = 0
        cov2_5_10 = 0
        cov2_1_5 = 0

        for key in target_amplicon_no_hp:
            if int(record_dir[key] - ntc_record_dir[key]) >= 50:
                # print(key, record_dir[key])
                cov2_50_plus += 1
            elif 20 <= int(record_dir[key] - ntc_record_dir[key]) < 50:
                cov2_20_50 += 1
            elif 10 <= int(record_dir[key] - ntc_record_dir[key]) < 20:
                cov2_10_20 += 1
            elif 5 <= int(record_dir[key] - ntc_record_dir[key]) < 10:
                cov2_5_10 += 1
            elif 1 <= int(record_dir[key] - ntc_record_dir[key]) < 5:
                cov2_1_5 += 1

        # print("cov2_50_plus", cov2_50_plus)
        # print("cov2_20_50", cov2_20_50)
        # print("cov2_10_20", cov2_10_20)
        # print("cov2_5_10", cov2_5_10)
        # print("cov2_1_5", cov2_1_5)

        if cov2_50_plus >= 3:
            result_mark = 'POS_3'
        elif cov2_50_plus == 2:
            result_mark = 'POS_2'
        elif cov2_50_plus == 1 and int(cov2_20_50 + cov2_10_20 + cov2_5_10) >= 2:
            result_mark = 'POS_2'
        elif int(cov2_50_plus + cov2_20_50) >= 1:
            result_mark = 'POS_1'
        elif int(cov2_50_plus + cov2_20_50 + cov2_10_20 + cov2_5_10) >= 2:
            result_mark = 'POS_1'
        elif int(cov2_50_plus + cov2_20_50 + cov2_10_20 + cov2_5_10 + cov2_1_5) >= 3:
            result_mark = 'POS_1'
        elif int(cov2_50_plus + cov2_20_50 + cov2_10_20 + cov2_5_10 + cov2_1_5) >= 1:
            # result_mark = 'POS_0'
            if int(record_dir[housekeeping] - ntc_record_dir[housekeeping]) >= 1000:
                result_mark = 'POS_0'
            else:
                result_mark = 'UNK'
        else:
            if int(record_dir[housekeeping] - ntc_record_dir[housekeeping]) >= 1000:
                result_mark = 'NEG'
            else:
                result_mark = 'UNK'
        # print(bc, '\t', result_mark)

        read_number = bc_df['read_number'].values[0]
        base_number = bc_df['total_base'].values[0]
        this_result = [bc, result_mark, read_number, base_number] + this_result_end
        this_result = '\t'.join([str(s) for s in this_result]) + '\n'
        summarized_result.append(this_result)

    with open(updated_result, 'a') as outfile:
        for line in summarized_result:
            outfile.writelines(line)
            logging.info(line.strip())
            # print(line.strip())


# summarize()
